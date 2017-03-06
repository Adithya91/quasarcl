#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <sstream>
#include <iostream>

#include <Rcpp.h>

namespace quasarcl {
	
	
		
	inline std::string readFile(std::string filename)
	{
		std::ifstream in;
		in.open(filename.c_str());
		if (in.is_open())
		{
			in.seekg(0, std::ios::end);
			size_t size = in.tellg();
			std::string buffer(size, ' ');
			in.seekg(0, std::ios::beg);
			in.read(&buffer[0], size);
			in.close();
			return buffer;
		}
		else
		{
			Rcpp::stop("Error occurred while opening file: " + filename);
		}
	}

	

	inline std::vector<std::string> getSources(std::vector<std::string>& filenames)
	{
		std::vector<std::string> sources;
		for (const auto& filename : filenames)
		{
			sources.push_back(readFile(filename));
		}
		return sources;
	}
	
	

	inline size_t calcGlobalSize(size_t maxWorkGroupSize, size_t dataSize)
	{
		size_t size = dataSize;
		size_t remainder = size % maxWorkGroupSize;
		if (remainder != 0)
		{
			size += maxWorkGroupSize - remainder;
		}
		if (size < dataSize)
		{
			Rcpp::stop("Error in calculating global_work_size.");
		}
		return size;
	}
}

#endif