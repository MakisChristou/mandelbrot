# cundelbrot: A CUDA based fractal renderer


## Compilation (Ubuntu 22.04)
```
sudo apt-get install nvidia-cuda-toolkit libsdl2-dev
cd Mandelbrot
make
```

## Running

```
./cundelbrot
```

## Hotkeys

| Key       | Function  | 
| ------------- |:-------------:|
| Arrow Keys | Pan around |
| W/S | +/- max iterations |
| E/D | +/- anti aliasing |
| Scrolwheel | Zoom In/Out |
| 0 | Reset to default |
| i | Set start keyframe |
| o | Set end keyframe |
| p | Playback animation |
| r | Render animation | 

## GPU vs CPU

| Image       | Ryzen 5950x (1 thread)  | RTX 3080 | RTX 3080 (kernel only)| 
| ------------- |:-------------:|:-------------:|:-------------:|
| Sample 1 | 0m12,828s | 0m2,275s | 88 ms  |
| Sample 2 | 2m19,883s | 0m3,274s | 1072 ms |
| Sample 3 | 0m30,843s | 0m2,596s | 290 ms | 
| Sample 4 | 1m30,134s | 0m3,363s | 1083 ms |


## Sample Images

### Sample1, Resolution: 5000x5000, n_max = 64, s_max = 8, range: -2.0, 2.0
![Sample 1](images/file1.png)

### Sample2, Resolution: 5000x5000, n_max = 256, s_max = 8, range: 0.2, 0.5
![Sample 2](images/file2.png)

### Sample3, Resolution: 5000x5000, n_max = 256, s_max = 8, range: 0.35, 0.45
![Sample 3](images/file3.png)

### Sample4, Resolution: 5000x5000, n_max = 512, s_max = 8, range: 0.35, 0.36
![Sample 3](images/file4.png)