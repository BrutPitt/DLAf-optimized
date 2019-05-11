#include <chrono>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

// number of DLA dimensions (must be 2 or 3)
#define DLA_DIM 3

// change to use float or double precision
using precisionT = float;


// Vector represents a point or a vector
/////////////////////////////////////////////////////////////////////
template <typename T> class Vector {
public:
    Vector<T>()              : m_X(0), m_Y(0), m_Z(0) {}
    Vector<T>(T x, T y)      : m_X(x), m_Y(y), m_Z(0) {}
    Vector<T>(T x, T y, T z) : m_X(x), m_Y(y), m_Z(z) {}

    T X() const { return m_X; }
    T Y() const { return m_Y; }
    T Z() const { return m_Z; }

    T Length() const {
        return std::sqrt(m_X * m_X + m_Y * m_Y + m_Z * m_Z);
    }

    T LengthSquared() const {
        return m_X * m_X + m_Y * m_Y + m_Z * m_Z;
    }

    T Distance(const Vector<T> &v) const {
        const T dx = m_X - v.m_X;
        const T dy = m_Y - v.m_Y;
        const T dz = m_Z - v.m_Z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    Vector<T> Normalized() const {
        const T m = 1 / Length();
        return Vector<T>(m_X * m, m_Y * m, m_Z * m);
    }

    Vector<T> operator+(const Vector<T> &v) const {
        return Vector<T>(m_X + v.m_X, m_Y + v.m_Y, m_Z + v.m_Z);
    }

    Vector<T> operator-(const Vector<T> &v) const {
        return Vector<T>(m_X - v.m_X, m_Y - v.m_Y, m_Z - v.m_Z);
    }

    Vector<T> operator*(const T a) const {
        return Vector<T>(m_X * a, m_Y * a, m_Z * a);
    }

    Vector<T> &operator+=(const Vector<T> &v) {
        m_X += v.m_X; m_Y += v.m_Y; m_Z += v.m_Z;
        return *this;
    }

private:
    T m_X, m_Y, m_Z;
};

// Pseudo Random Generator to use
/////////////////////////////////////////////////////////////////////
#define DLAF_USE_FAST_RANDOM    //comment to use std::mt19937 & std::uniform_real_distribution
#ifdef DLAF_USE_FAST_RANDOM
//Use Marsaglia fast random generator
    #include "fastRandom.h"
    using namespace fstRnd;
    //#define DLAF_USE_64BIT_GENERATOR
    #ifdef DLAF_USE_64BIT_GENERATOR
        floatfastRandomClass<precisionT, fastRand64> fastRandom;    
    #else
        floatfastRandomClass<precisionT, fastRand32> fastRandom;
    #endif    
    #define DLA_RANDOM_NORM fastRandom.VNI()    // [-1.0, 1.0]
    #define DLA_RANDOM_01   fastRandom.UNI()    // [ 0.0, 1.0]
#else
// std::mt19937 & std::uniform_real_distribution
// Random returns a uniformly distributed random number between lo and hi
    #include <random>
    template<typename T> T Random(const T lo = 0, const T hi = 1) {
        static thread_local std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<T> dist(lo, hi);
        return dist(gen);
    }
    #define DLA_RANDOM_NORM Random(T(-1.0), T(1.0))
    #define DLA_RANDOM_01   Random(T( 0.0), T(1.0))
#endif


// Library to use for spatial index
/////////////////////////////////////////////////////////////////////
#define DLAF_USE_FLANN_LIBRARY   // comment to use boost instead nanoflann
#ifdef DLAF_USE_FLANN_LIBRARY
// nanoflann is used for its spatial index
    #include "nanoflann.hpp"

    #define parentPOINT(PARENT) m_Points.pts[PARENT]
    #define thisPOINT m_Points.pts

template <typename T> struct pointCloud
{
    std::vector<Vector<T>> pts;
    
    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }
    
    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const
    { return dim==0 ? pts[idx].X() : (dim==1 ? pts[idx].Y() : pts[idx].Z()); }
    
    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX> bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
    
};

using tPointCloud =  pointCloud<precisionT>;
using tKDTreeDistanceFunc = nanoflann::L2_Adaptor<precisionT, tPointCloud>;
using tKDTree = nanoflann::KDTreeSingleIndexDynamicAdaptor<tKDTreeDistanceFunc, tPointCloud, DLA_DIM>;
#else
// boost is used for its spatial index
    #include <boost/function_output_iterator.hpp>
    #include <boost/geometry/geometry.hpp>

    using BoostPoint = boost::geometry::model::point<precisionT, DLA_DIM, boost::geometry::cs::cartesian>;
    using IndexValue = std::pair<BoostPoint, uint32_t>;
    using Index = boost::geometry::index::rtree<IndexValue, boost::geometry::index::linear<4>>;

    #define parentPOINT(PARENT) m_Points[PARENT]
    #define thisPOINT m_Points
#endif

// Lerp linearly interpolates from a to b by distance.
template<typename T> Vector<T> Lerp(const Vector<T> &a, const Vector<T> &b, const T d) {
    return a + (b - a).Normalized() * d;
}

// default parameters (documented below)
const precisionT DefaultParticleSpacing = 1;
const precisionT DefaultAttractionDistance = 3;
const precisionT DefaultMinMoveDistance = 1;
const precisionT DefaultStickiness = 1;
const int DefaultStubbornness = 0;

// Model holds all of the particles and defines their behavior.
template<typename T> class Model {
public:
    Model<T>() : 
        m_ParticleSpacing(DefaultParticleSpacing),
        m_AttractionDistance(DefaultAttractionDistance),
        m_MinMoveDistance(DefaultMinMoveDistance),
        m_Stubbornness(DefaultStubbornness),
        m_Stickiness(DefaultStickiness),
        m_BoundingRadius(0)        
#ifdef  DLAF_USE_FLANN_LIBRARY
        { m_Index = new tKDTree(DLA_DIM, m_Points, nanoflann::KDTreeSingleIndexAdaptorParams(10 /*max leaf*/)); }
        ~Model<T>() { delete m_Index; }
#else        
        {}
#endif        

    void SetParticleSpacing(const T a) {
        m_ParticleSpacing = a;
    }

    void SetAttractionDistance(const T a) {
        m_AttractionDistance = a;
    }

    void SetMinMoveDistance(const T a) {
        m_MinMoveDistance = a;
    }

    void SetStubbornness(const int a) {
        m_Stubbornness = a;
    }

    void SetStickiness(const T a) {
        m_Stickiness = a;
    }

#ifdef DLAF_USE_FLANN_LIBRARY
    // Add adds a new particle with the specified parent particle
    void Add(const Vector<T>& p, const int parent = -1) {
        size_t id = m_Points.pts.size();
        m_Points.pts.push_back(p);
        m_JoinAttempts.push_back(0);
        m_Index->addPoints(id, id);

        m_BoundingRadius = std::max(m_BoundingRadius, p.Length() + m_AttractionDistance);
        std::cout << id << "," << parent << "," << p.X() << "," << p.Y() << "," << p.Z() << std::endl;

    }

    // Nearest returns the index of the particle nearest the specified point
    uint32_t Nearest(const Vector<T> &point) const {
        size_t ret_index;
        T out_dist_sqr = m_AttractionDistance;
        nanoflann::KNNResultSet<T> resultSet(1);
        resultSet.init(&ret_index, &out_dist_sqr );
        m_Index->findNeighbors(resultSet, (const T *) &point, 
                               nanoflann::SearchParams(0, //how many leafs to visit - not used in nanoflann
                                                       m_AttractionDistance*.5)); //search for eps-approximate neighbours
        return ret_index;
    }
#else
    // Add adds a new particle with the specified parent particle
    void Add(const Vector<T> &p, const int parent = -1) {
        const uint32_t id = m_Points.size();
        m_Index.insert(std::make_pair(BoostPoint(p.X(), p.Y(), p.Z()), id));
        m_Points.push_back(p);
        m_JoinAttempts.push_back(0);
        m_BoundingRadius = std::max(m_BoundingRadius, p.Length() + m_AttractionDistance);
        std::cout << id << "," << parent << "," << p.X() << "," << p.Y() << "," << p.Z() << std::endl;
    }
    // Nearest returns the index of the particle nearest the specified point
    uint32_t Nearest(const Vector<T> &p) const {
        uint32_t result = -1;
        m_Index.query(
            boost::geometry::index::nearest(BoostPoint(p.X(), p.Y(), p.Z()), 1),
            boost::make_function_output_iterator([&result](const auto &value) {
                result = value.second;
            }));
        return result;
    }
#endif
    // RandomInUnitSphere returns a random, uniformly distributed point inside the
    // unit sphere (radius = 1)
    Vector<T> RandomInUnitSphere() const {
        Vector<T> p;
        do {
            p = Vector<T>(DLA_RANDOM_NORM, 
                          DLA_RANDOM_NORM, 
                          DLA_DIM == 2 ? T(0) : DLA_RANDOM_NORM);
        } while(p.Length() >= T(1.0));

        return p;
    }

    // RandomStartingPosition returns a random point to start a new particle
    Vector<T> RandomStartingPosition() const {
        const T d = m_BoundingRadius;
        return RandomInUnitSphere().Normalized() * d;
    }

    // ShouldReset returns true if the particle has gone too far away and
    // should be reset to a new random starting position
    bool ShouldReset(const Vector<T> &p) const {
        return p.Length() > m_BoundingRadius * T(2);
    }

    // ShouldJoin returns true if the point should attach to the specified
    // parent particle. This is only called when the point is already within
    // the required attraction distance.
    bool ShouldJoin(const Vector<T> &p, const int parent) {
        return (m_JoinAttempts[parent]++ < m_Stubbornness) ? false : DLA_RANDOM_01 <= m_Stickiness;
    }

    // PlaceParticle computes the final placement of the particle.
    Vector<T> PlaceParticle(const Vector<T> &p, const int parent) const {
        return Lerp(parentPOINT(parent), p, m_ParticleSpacing);
    }

    // MotionVector returns a vector specifying the direction that the
    // particle should move for one iteration. The distance that it will move
    // is determined by the algorithm.
    Vector<T> MotionVector(const Vector<T> &p) const {
        return RandomInUnitSphere();
    }

    // AddParticle diffuses one new particle and adds it to the model
    Vector<T> AddParticle() {
        // compute particle starting location
        Vector<T> p = RandomStartingPosition();

        // do the random walk
        while (true) {
            // get distance to nearest other particle
            const int parent = Nearest(p);
            const T d = p.Distance(parentPOINT(parent));

            // check if close enough to join
            if (d < m_AttractionDistance) {
                if (!ShouldJoin(p, parent)) {
                    // push particle away a bit
                    p = Lerp(parentPOINT(parent), p, m_AttractionDistance + m_MinMoveDistance);
                    continue;
                }

                // adjust particle position in relation to its parent
                p = PlaceParticle(p, parent);

                // adjust particle pos in relation to its parent and add the point
                Add(PlaceParticle(p, parent), parent);
                return thisPOINT.back();
            }

            // move randomly
            const T m = std::max(m_MinMoveDistance, d - m_AttractionDistance);
                p += MotionVector(p).Normalized() * m;

            // check if particle is too far away, reset if so
            if (ShouldReset(p)) p = RandomStartingPosition();
        }
    }

private:
    // m_ParticleSpacing defines the distance between particles that are
    // joined together
    T m_ParticleSpacing;

    // m_AttractionDistance defines how close together particles must be in
    // order to join together
    T m_AttractionDistance;

    // m_MinMoveDistance defines the minimum distance that a particle will move
    // during its random walk
    T m_MinMoveDistance;

    // m_Stubbornness defines how many interactions must occur before a
    // particle will allow another particle to join to it.
    int m_Stubbornness;

    // m_Stickiness defines the probability that a particle will allow another
    // particle to join to it.
    T m_Stickiness;

    // m_BoundingRadius defines the radius of the bounding sphere that bounds
    // all of the particles
    T m_BoundingRadius;

    // m_JoinAttempts tracks how many times other particles have attempted to
    // join with each finalized particle
    std::vector<int> m_JoinAttempts;


#ifdef DLAF_USE_FLANN_LIBRARY
    // m_Index is the spatial index used to accelerate nearest neighbor queries
    tKDTree *m_Index;
    // m_Points stores the final particle positions
    tPointCloud m_Points;
#else
    // m_Index is the spatial index used to accelerate nearest neighbor queries
    Index m_Index;
    // m_Points stores the final particle positions
    std::vector<Vector<T>> m_Points;
#endif
};

int main() {
    // use float or double precision
    //  using precisionT = float;

    // create the model
    Model<precisionT> model;

    // add seed point(s)
    model.Add(Vector<precisionT>());

    // {
    //     const int n = 3600;
    //     const precisionT r = 1000;
    //     for (int i = 0; i < n; i++) {
    //         const precisionT t = (precisionT)i / n;
    //         const precisionT a = t * 2 * M_PI;
    //         const precisionT x = std::cos(a) * r;
    //         const precisionT y = std::sin(a) * r;
    //         model.Add(Vector<precisionT>(x, y, 0));
    //     }
    // }

    auto start = std::chrono::high_resolution_clock::now();
    // run diffusion-limited aggregation
    for (int i = 0; i < 1000000; i++) {
        model.AddParticle();
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;

    std::cerr << "Time elapsed: " << diff.count() << " sec." << std::endl;

    return 0;
}
