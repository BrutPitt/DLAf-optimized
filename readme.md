# DLAf-optimized

This is an optimized version of [Michael Fogleman's dlaf](https://github.com/fogleman/dlaf) project that I used in [**glChAoS.P / wglChAoS.P**](https://github.com/BrutPitt/glChAoS.P) for large-scale DLA3D growth.

[![](https://raw.githubusercontent.com/BrutPitt/myRepos/master/glChAoSP/screenShots/dla3D.jpg)](https://youtu.be/_R7akPh64XU)
Video Real-time rendering (80 FPS) of 6M points/particles: https://youtu.be/_R7akPh64XU 

All the changes has been made in [**glChAoS.P / wglChAoS.P**](https://github.com/BrutPitt/glChAoS.P), but I thought I'd make them available here, in a separate repository, for simplicity.

### Main changes
- The use of *single file header* [**nanoflann**](https://github.com/jlblancoc/nanoflann) (Fast Library for Approximate Nearest Neighbors) for spatial index, instead of boost library.
- The use of [**fastPRNG**](https://github.com/BrutPitt/fastPRNG), based xoshiro / xoroshiro (Blackman/Vigna), xorshift and other Marsaglia algorithms.
- The use of `template class` for `Vector` and `Model` with possibilty to use single or double precison floating points.

The [Michael Fogleman](https://github.com/fogleman/dlaf)'s original option (such as boost library, or std::mt19937 random generator) are anyway aviables via internal `#define` 

**Performance and differences calculating 1M points**

Tests were performed on: Intel I5 6500 - 3.2 GHz - single thread - DDR4 2666MHz - source code compiled with gcc 8.3.1 on Linux Fedora29.

- **Single precision tests** (float)
  - nanoflann / fastRandom:   **63.21 sec.**
  - nanoflann / std::mt19937: **65.04 sec.**
  - boost lib / fastRandom: **123.78 sec.**
  - boost lib / std::mt19937: **125.33 sec.**
- **Double precision tests** (double)
  - nanoflann / fastRandom:   **68.14 sec.**
  - nanoflann / std::mt19937: **70.06 sec.**
  - boost lib / fastRandom: **124.62 sec.**
  - boost lib / std::mt19937: **126.57 sec.**

### External Dependences

**NONE**: all files are already included in the project/repository

- `nanoflann.h` -> header-only library for KD-Trees of datasets point clouds: [**nanoflann**](https://github.com/jlblancoc/nanoflann)
- `fastPRNG.h` -> fast PRNG (Pseudo-Random Number Generator) based xoshiro / xoroshiro (Blackman/Vigna), xorshift and other Marsaglia algorithms.: [**fastPRNG**](https://github.com/BrutPitt/fastPRNG)

If you want use Boost library, like the original project, it's not included in the repository, but can be downloaded from [https://www.boost.org/](https://www.boost.org/)
It's not necessary to build the library, only headers files is enough. (about 160MB) 

### Output Format
This program produces only a simple text output, for both for 2D and 3D, always in the form:
```bash
0,-1,0,0,0
1,0,-0.191125,-0.446429,0.87417
2,1,0.697415,-0.821024,1.13908
3,0,-0.464858,0.327307,-0.822665
4,2,1.63213,-0.595306,0.86457
5,1,-1.06683,-0.335015,1.34398
6,5,-1.88992,-0.368307,1.91091
7,6,-2.20384,0.449373,2.39347
...
```
Where the columns are: point_id, parent_id, x, y, z (for 2D z is always 0)

You can redirect the output in a file with command: `./dlaf > output.csv`

If you instead want to export a DLA3D object in PLY file format (ascii/binary), please use [**glChAoS.P / wglChAoS.P**](https://github.com/BrutPitt/glChAoS.P); you can also import a PLY object, where to grow the DLA:

| [**glChAoS.P**](https://github.com/BrutPitt/glChAoS.P) - Thomas Attractor |  [**glChAoS.P**](https://github.com/BrutPitt/glChAoS.P) - DLA on Thomas Attractor |
| :-----: | :----: |
|[![](https://user-images.githubusercontent.com/16171743/57495764-ebb7d180-72cf-11e9-9da6-b2f4be227951.jpg)](https://user-images.githubusercontent.com/16171743/57495764-ebb7d180-72cf-11e9-9da6-b2f4be227951.jpg) | [![](https://user-images.githubusercontent.com/16171743/57906343-e6cbc280-7879-11e9-8ada-fc639dce60d4.jpg)](https://user-images.githubusercontent.com/16171743/57906343-e6cbc280-7879-11e9-8ada-fc639dce60d4.jpg)|

Watch **[DLA3D over Thomas Attractor Video](https://youtu.be/Eu69KcPms_s)**



### Build DLAf-optimized

Just clone it form github or download the archive and decompress it in a folder whatever and  use `make` from it.
Obviously need to install the GCC C/C++ v.6 or higher (or clang v.6 or higher) compilers and associated tools (such as *make*), or other C++14 compiler.

To install gcc C/C++ in Linux:
* Debian, Ubuntu: `sudo apt install build-essential`
* Fedora, RedHat: `sudo dnf install make gcc-c++`

In the `dlaf` folder there is a very simple `cmake` file from wich you can build also a *project solution* for Microsoft Visual Studio.


### DLAf-optimized - Hooks & Parameters 
Description form [Michael Fogleman](https://github.com/fogleman/dlaf)'s page:

The code implements a standard diffusion-limited aggregation algorithm. But there are several parameters and code hooks that let you tweak its behavior.

The following parameters can be set on a `Model` instance.

| Parameter | Description |
| --- | --- |
| `AttractionDistance` | Defines how close together particles must be in order to join together. |
| `ParticleSpacing` | Defines the distance between particles when they become joined together. |
| `MinMoveDistance` | Defines the minimum distance that a particle will move in an iteration during its random walk. |
| `Stubbornness` | Defines how many join attempts must occur before a particle will allow another particle to join to it. |
| `Stickiness` | Defines the probability that a particle will allow another particle to join to it. |

The following hooks allow you to define the algorithm behavior in small, well-defined functions.

| Hook | Description |
| --- | --- |
| `RandomStartingPosition()` | Returns a starting position for a new particle to begin its random walk. |
| `ShouldReset(p)` | Returns true if the particle has gone too far away and should be reset to a new random starting position. |
| `ShouldJoin(p, parent)` | Returns true if the point should attach to the specified parent particle. This is only called when the point is already within the required attraction distance. If false is returned, the particle will continue its random walk instead of joining to the other particle. |
| `PlaceParticle(p, parent)` | Returns the final placement position of the particle. |
| `MotionVector(p)` | Returns a vector specifying the direction that the particle should move for one iteration. The distance that it will move is determined by the algorithm. |



