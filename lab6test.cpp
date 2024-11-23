#include "lab6.h"
#include <chrono>
#include <fstream>
#include "C:\Users\Влад\source\repos\Cubic interpolation spline\Cubic interpolation spline.cpp"

int main() {

	Fourier Obj = Fourier();

	std::vector<std::complex<double>> signal_dft;

	for (int i = 0; i < 512; i++) {

		signal_dft.push_back(Obj.generate(i));

	}

	std::vector<std::complex<double>> signal_fft = signal_dft;

	std::cout << "\n\nSIGNAL:\n\n";

	for (int i = 0; i < 512; i++) {

		std::cout << i << ' ' << signal_dft.at(i) << std::endl;

	}

	auto s = std::chrono::high_resolution_clock::now();

	Obj.DFT(signal_dft);

	auto e = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time = e - s;
	std::cout << "\n\nExecution time (DFT): " << elapsed_time.count() << " seconds\n\n";

	for (int i = 0; i < 512; i++) {

		if (Obj.amplitude(signal_dft.at(i)) > 1e-9) std::cout << i << ' ' << signal_dft.at(i) << ' ' << Obj.amplitude(signal_dft.at(i)) << ' ' << Obj.phase(signal_dft.at(i)) << std::endl;

	}

	s = std::chrono::high_resolution_clock::now();

	Obj.IDFT(signal_dft);

	e = std::chrono::high_resolution_clock::now();
	elapsed_time = e - s;
	std::cout << "\n\nExecution time (IDFT): " << elapsed_time.count() << " seconds\n\n";

	for (int i = 0; i < 512; i++) {

		std::cout << i << ' ' << signal_dft.at(i) << std::endl;

	}

	//FFT

	s = std::chrono::high_resolution_clock::now();

	Obj.FFT(signal_fft);

	e = std::chrono::high_resolution_clock::now();
	elapsed_time = e - s;
	std::cout << "\n\nExecution time (FFT): " << elapsed_time.count() << " seconds\n\n";

	for (int i = 0; i < 512; i++) {

		if (Obj.amplitude(signal_fft.at(i)) > 1e-9) std::cout << i << ' ' << signal_fft.at(i) << ' ' << Obj.amplitude(signal_fft.at(i)) << ' ' << Obj.phase(signal_fft.at(i)) << std::endl;

	}

	s = std::chrono::high_resolution_clock::now();

	Obj.IDFT(signal_fft);

	e = std::chrono::high_resolution_clock::now();
	elapsed_time = e - s;
	std::cout << "\n\nExecution time (IFFT): " << elapsed_time.count() << " seconds\n\n";

	for (int i = 0; i < 512; i++) {

		std::cout << i << ' ' << signal_fft.at(i) << std::endl;

	}

	//генерация шума

	std::cout << "\n\nWITH NOISE:\n\n";

	std::vector<std::complex<double>> signal_noise;

	for (int i = 0; i < 512; i++) {

		signal_noise.push_back(Obj.generate_noise(i));

		std::cout << i << ' ' << signal_noise.back() << std::endl;

	}

	//сплайн

	std::vector<std::pair<double, double>> WITH_GRID;

	for (int i = 0; i < 512; i++) {

		WITH_GRID.push_back(std::pair<double, double>(static_cast<double>(i), signal_noise.at(i).real()));

	}

	Spline WITH(1., 0., 1., 512, WITH_GRID);

	WITH.gMake();

	std::ofstream ofs("withnoise.txt");

	for (int i = 0; i < 512; i++) {

		ofs << static_cast<double>(i) << ' ' << WITH.g(static_cast<double>(i)) << std::endl;

		ofs << static_cast<double>(i) + 0.25 << ' ' << WITH.g(static_cast<double>(i) + 0.25) << std::endl;

		ofs << static_cast<double>(i) + 0.5 << ' ' << WITH.g(static_cast<double>(i) + 0.5) << std::endl;

		ofs << static_cast<double>(i) + 0.75 << ' ' << WITH.g(static_cast<double>(i) + 0.75) << std::endl;

	}

	ofs.close(); ofs.clear();


	//DFT, обнуление, IDFT

	std::cout << "\n\nDFT WITH NOISE:\n\n";

	Obj.DFT(signal_noise);

	for (int i = 0; i < 512; i++) {

		if (Obj.amplitude(signal_noise.at(i)) > 1e-9) std::cout << i << ' ' << signal_noise.at(i) << ' ' << Obj.amplitude(signal_noise.at(i)) << ' ' << Obj.phase(signal_noise.at(i)) << std::endl;

	}

	signal_noise.at(229) = (0., 0.); signal_noise.at(283) = (0., 0.);

	std::cout << "\n\nIDFT WITHOUT NOISE:\n\n";

	Obj.IDFT(signal_noise);

	for (int i = 0; i < 512; i++) {

		std::cout << i << ' ' << signal_noise.at(i) << std::endl;

	}

	//сплайн

	std::vector<std::pair<double, double>> WITHOUT_GRID;

	for (int i = 0; i < 512; i++) {

		WITHOUT_GRID.push_back(std::pair<double, double>(static_cast<double>(i), signal_noise.at(i).real()));

	}

	Spline WITHOUT(1., 0., 1., 512, WITHOUT_GRID);

	WITHOUT.gMake();

	ofs.open("withoutnoise.txt");

	for (int i = 0; i < 512; i++) {

		ofs << static_cast<double>(i) << ' ' << WITHOUT.g(static_cast<double>(i)) << std::endl;

		ofs << static_cast<double>(i) + 0.25 << ' ' << WITHOUT.g(static_cast<double>(i) + 0.25) << std::endl;

		ofs << static_cast<double>(i) + 0.5 << ' ' << WITHOUT.g(static_cast<double>(i) + 0.5) << std::endl;

		ofs << static_cast<double>(i) + 0.75 << ' ' << WITHOUT.g(static_cast<double>(i) + 0.75) << std::endl;

	}

	ofs.close(); ofs.clear();

	return 0;

}