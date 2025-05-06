# simple-pgo
C code Library for studying pose graph optimization

## Dependency
```
sudo apt install libglfw3-dev libgl1-mesa-dev libglu1-mesa-dev libglew-dev pkg-config
```
## Build & Run
```
mkdir build && cd build
cmake ..
make
./simple_pgo
```

## File Tree
```
simple-pgo/
├── CMakeLists.txt
├── core/
│   ├── graph.c
│   ├── graph.h
│   ├── optimizer.c
│   ├── optimizer.h
│   ├── jacobian_ops.c
│   ├── jacobian_ops.h
│   └── CMakeLists.txt
├── types/
│   ├── vertex_se2.c
│   ├── vertex_se2.h
│   ├── edge_se2.c
│   ├── edge_se2.h
│   ├── vertex_se3.c
│   ├── edge_se3.c
│   └── CMakeLists.txt
├── include/
│   └── simple_pgo/
│       └── core_api.h
├── visualization/
│   ├── viewer.c              # OpenGL based visualization code
│   ├── viewer.h              # OpenGL based visualization code
│   ├── g2o_parser.c
│   ├── g2o_parser.h
│   └── CMakeLists.txt
├── main.c                    # execution file
└── data/                     # test data
```
