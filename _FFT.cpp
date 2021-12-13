#include <complex>
#include <vector>
#include <math.h>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

using namespace std;

namespace py = pybind11;

const double pi = acos(-1.0);

unsigned bit_reverse (unsigned a, int len) {
    a = ((a&0x55555555U)<<1)  | ((a&0xAAAAAAAAU)>>1);
    a = ((a&0x33333333U)<<2)  | ((a&0xCCCCCCCCU)>>2);
    a = ((a&0x0F0F0F0FU)<<4)  | ((a&0xF0F0F0F0U)>>4);
    a = ((a&0x00FF00FFU)<<8)  | ((a&0xFF00FF00U)>>8);
    a = ((a&0x0000FFFFU)<<16) | ((a&0xFFFF0000U)>>16);
    return a >> (32-len);
}

vector<complex<double>> fft(vector<complex<double>> in) {
	int in_len = in.size();
	int N = 1;
	while (N < in_len) {
		N <<= 1;
	}
	vector<complex<double>> out(N);
	for (int i = in_len; i < N; i++) {
		in.push_back(complex<double>(0, 0));
	}
	int bitlen = log2(N);
	for (int i = 0; i < N; ++i) {
		out.at(bit_reverse(i, bitlen)) = in.at(i);
	}
	for (int step = 2, mh = 1; step <= N; step <<= 1, mh <<= 1) {
		for (int i = 0; i < mh; ++i) {
			complex<double> wi = exp(complex<double>(0, i * pi / mh));
			for (int j = i, k = i + mh; j < N; j += step, k += step) {
				complex<double> u = out.at(j);
				complex<double> t = wi * out.at(k);
				out.at(j) = u + t;
				out.at(k) = u - t;
			}
		}
	}
	out.resize(in_len);
	return out;
}


vector<complex<double>> ifft(vector<complex<double>> in) {
	int in_len = in.size();
	int N = 1;
	while (N < in_len) {
		N <<= 1;
	}
	for (int i = in_len; i < N; i++) {
		in.push_back(complex<double>(0, 0));
	}
	vector<complex<double>> out(N);
	int bitlen = log2(N);
	for (int i = 0; i < N; ++i) {
		out.at(bit_reverse(i, bitlen)) = in.at(i);
	}
	for (int step = 2, mh = 1; step <= N; step <<= 1, mh <<= 1) {
		for (int i = 0; i < mh; ++i) {
			complex<double> wi = exp(complex<double>(0, -1 * i * pi / mh));
			for (int j = i, k = i + mh; j < N; j += step, k += step) {
				complex<double> u = out.at(j);
				complex<double> t = wi * out.at(k);
				out.at(j) = u + t;
				out.at(k) = u - t;
			}
		}
	}
	for (int i = 0; i < N; ++i) {
		out.at(i) /= N;
	}
	out.resize(in_len);
	return out;
}

PYBIND11_MODULE(_FFT, m) {
  m.doc() = "FFT";
  m.def("fft", &fft);
  m.def("ifft", &ifft);
}