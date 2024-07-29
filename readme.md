# HyperCutter

This repository contains the personel implementation of 

> **Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs**<br/>
> Julius Nehring-Wirxel, Philip Trettner, Leif Kobbelt <br/>
> [PDF](https://arxiv.org/pdf/2103.02486)

## Screenshot
### random plane sphere cut
![randomplanecut](screenshot/randomsphere00.png)
### mesh 2 bsp convertion
![mesh2bsp](screenshot/mesh2bsp00.png)
### virtual milling example

<video src="https://raw.githubusercontent.com/quadmotor/hypercut/main/screenshot/bunny.mp4"></video>
![bunny.mp4](screenshot/bunny.gif)

## Compile
```shell
git clone --recursive git@github.com:quadmotor/hypercut.git
cd hypercut
mkdir build
cd build
cmake ..
cmake --build .
```

## Usage

### random plane sphere cut
```shell
hypercuttertest.exe --randomsphere   --output  C:\tmp\hc\randoms
```
### mesh 2 bsp convertion
```shell
hypercutter.exe --mesh2bsp --input drill-trace1_1.hh.sat.obj   --output  C:\tmp\hc\randoms.obj
```
### vertical mill example
```shell
hypercutter.exe --mill --millresolution 30 --input bunny2.obj   --output  C:\tmp\hc\randoms.obj
```


## Limitations
This repo only implement the bsp tree and limit to 128bit.

## Credit
* meshlab
* https://polyscope.run/py/ 
* https://github.com/pmp-library/pmp-library
* https://github.com/devatrun/slimcpplib
* https://github.com/wjakob/instant-meshes/
* AABB-triangle overlap test code  

