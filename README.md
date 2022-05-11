# Mandelbrot set renderer


## Compilation

```
cd Mandelbrot
make
```

## Running

```
./mandelbrot > file.ppm
```

## Converting output to png

```
sudo apt install sudo apt install imagemagick
convert file.ppm file.png
```

## Sample Images

### Resolution: 1000x1000, n_max = 128, s_max = 1, range: -2.0, 2.0
![Sample 1](images/file1.png)

### Resolution: 1000x1000, n_max = 256, s_max = 1, range: 0.2, 0.5
![Sample 2](images/file2.png)

### Resolution: 1000x1000, n_max = 256, s_max = 1, range: 0.35, 0.45
![Sample 3](images/file3.png)

### Resolution: 1000x1000, n_max = 512, s_max = 1, range: 0.35, 0.36
![Sample 3](images/file4.png)