#### Linux

```bash
g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` _FFT.cpp -o _FFT`python3-config --extension-suffix`
```

#### Mac

```bash
g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` _FFT.cpp -o _FFT`python3-config --extension-suffix` 
```

