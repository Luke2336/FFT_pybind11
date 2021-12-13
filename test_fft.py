import _FFT

a = [1, 2, 0, 0]
b = [2, 3, 0, 0]
a_fft = _FFT.fft(a)
b_fft = _FFT.fft(b)
print(a_fft, b_fft)
ab_fft = []
for i in range(len(a_fft)):
	ab_fft.append(a_fft[i] * b_fft[i])
ab = _FFT.ifft(ab_fft)
print(ab)