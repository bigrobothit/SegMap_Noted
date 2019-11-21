// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2012,
François Pomerleau and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <f dot pomerleau at gmail dot com> and
<stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "IO.h"
#include "IOFunctions.h"
#include "InspectorsImpl.h"

// For logging
#include "PointMatcherPrivate.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ctype.h>
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/foreach.hpp"

#ifdef WIN32
#define strtok_r strtok_s
#endif // WIN32

namespace PointMatcherSupport {
	namespace {
		const int one = 1;
	}
	const bool isBigEndian = *reinterpret_cast<const unsigned char*>(&one) == static_cast<unsigned char>(0);
	const int oneBigEndian = isBigEndian ? 1 : 1 << 8 * (sizeof(int) - 1);
}

using namespace std;
using namespace PointMatcherSupport;


// Tokenize a string, excepted if it begins with a '%' (a comment in CSV)
static std::vector<string> csvLineToVector(const char* line)
{
	std::vector<string> parsedLine;
	char delimiters[] = " \t,;";
	char *token;
	char tmpLine[1024];
	char *brkt = 0;
	strcpy(tmpLine, line);
	token = strtok_r(tmpLine, delimiters, &brkt);
	if(line[0] != '%') // Jump line if it's commented
	{
		while (token)
		{
			parsedLine.push_back(string(token));
			token = strtok_r(NULL, delimiters, &brkt);
		}
	}

	return parsedLine;
}

// Open and parse a CSV file, return the data
CsvElements parseCsvWithHeader(const std::string& fileName)
{
	validateFile(fileName);
	
	ifstream is(fileName.c_str());

	unsigned elementCount=0;
	std::map<string, unsigned> keywordCols;
	CsvElements data;

	bool firstLine(true);
	unsigned lineCount=0;
	string line;
	//while (!is.eof())
	while (std::getline(is, line))
	{
		//char line[1024];
		//is.getline(line, sizeof(line));
		//line[sizeof(line)-1] = 0;



		if(firstLine)
		{
			std::vector<string> header = csvLineToVector(line.c_str());
				
			elementCount = header.size();
			for(unsigned int i = 0; i < elementCount; i++)
			{
				keywordCols[header[i]] = i;
			}

			firstLine = false;
		}
		else // load the rest of the file
		{
			std::vector<string> parsedLine = csvLineToVector(line.c_str());
			if(parsedLine.size() != elementCount && parsedLine.size() !=0)
			{
				stringstream errorMsg;
				errorMsg << "Error at line " << lineCount+1 << ": expecting " << elementCount << " columns but read " << parsedLine.size() << " elements.";
				throw runtime_error(errorMsg.str());	
			}

			for(unsigned int i = 0; i < parsedLine.size(); i++)
			{
				for(BOOST_AUTO(it,keywordCols.begin()); it!=keywordCols.end(); it++)
				{
					if(i == (*it).second)
					{
						data[(*it).first].push_back(parsedLine[i]);	
					}
				}
			}
		}

		lineCount++;
	}
	
	// Use for debug
	
	//for(BOOST_AUTO(it,data.begin()); it!=data.end(); it++)
	//{
	//	cout << "--------------------------" << endl;
	//	cout << "Header: |" << (*it).first << "|" << endl;
	//	//for(unsigned i=0; i<(*it).second.size(); i++)
	//	//{
	//	//	cout << (*it).second[i] << endl;
	//	//}
	//}
	

	return data;
}


//! Constructor, leave fields blank if unused
template<typename T>
PointMatcherIO<T>::FileInfo::FileInfo(const std::string& readingFileName, const std::string& referenceFileName, const std::string& configFileName, const TransformationParameters& initialTransformation, const TransformationParameters& groundTruthTransformation, const Vector& gravity):
	readingFileName(readingFileName),
	referenceFileName(referenceFileName),
	configFileName(configFileName),
	initialTransformation(initialTransformation),
	groundTruthTransformation(groundTruthTransformation),
	gravity(gravity)
{}

template struct PointMatcherIO<float>::FileInfo;
template struct PointMatcherIO<double>::FileInfo;

// Empty constructor
template<typename T>
PointMatcherIO<T>::FileInfoVector::FileInfoVector()
{
}

//! Load a vector of FileInfo from a CSV file.
/**
	@param fileName name of the CSV file
	@param dataPath path relative to which the point cloud CSV or VTK will be resolved
	@param configPath path relative to which the yaml configuration files will be resolved
	
	The first line of the CSV file must contain a header. The supported tags are:
	- reading: file name of the reading point cloud
	- reference: file name of the reference point cloud
	- config: file name of the YAML configuration of the ICP chain
	- iTxy: initial transformation, coordinate x,y
	- gTxy: ground-truth transformation, coordinate x,y
	Note that the header must at least contain "reading".
*/
template<typename T>
PointMatcherIO<T>::FileInfoVector::FileInfoVector(const std::string& fileName, std::string dataPath, std::string configPath)
{
	if (dataPath.empty())
	{
		#if BOOST_FILESYSTEM_VERSION >= 3
		dataPath = boost::filesystem::path(fileName).parent_path().string();
		#else
		dataPath = boost::filesystem::path(fileName).parent_path().file_string();
		#endif
	}
	if (configPath.empty())
	{
		#if BOOST_FILESYSTEM_VERSION >= 3
		configPath = boost::filesystem::path(fileName).parent_path().string();
		#else
		configPath = boost::filesystem::path(fileName).parent_path().file_string();
		#endif
	}
	
	const CsvElements data = parseCsvWithHeader(fileName);
	
	// Look for transformations
	const bool found3dInitialTrans(findTransform(data, "iT", 3));
	bool found2dInitialTrans(findTransform(data, "iT", 2));
	const bool found3dGroundTruthTrans(findTransform(data, "gT", 3));
	bool found2dGroundTruthTrans(findTransform(data, "gT", 2));
	if (found3dInitialTrans)
		found2dInitialTrans = false;
	if (found3dGroundTruthTrans)
		found2dGroundTruthTrans = false;
	
	// Check for consistency
	if (found3dInitialTrans && found2dGroundTruthTrans)
		throw runtime_error("Initial transformation is in 3D but ground-truth is in 2D");
	if (found2dInitialTrans && found3dGroundTruthTrans)
		throw runtime_error("Initial transformation is in 2D but ground-truth is in 3D");
	CsvElements::const_iterator readingIt(data.find("reading"));
	if (readingIt == data.end())
		throw runtime_error("Error transfering CSV to structure: The header should at least contain \"reading\".");
	CsvElements::const_iterator referenceIt(data.find("reference"));
	CsvElements::const_iterator configIt(data.find("config"));
	
	// Load reading
	const std::vector<string>& readingFileNames = readingIt->second;
	const unsigned lineCount = readingFileNames.size();
	boost::optional<std::vector<string> > referenceFileNames;
	boost::optional<std::vector<string> > configFileNames;
	if (referenceIt != data.end())
	{
		referenceFileNames = referenceIt->second;
		assert (referenceFileNames->size() == lineCount);
	}
	if (configIt != data.end())
	{
		configFileNames = configIt->second;
		assert (configFileNames->size() == lineCount);
	}

	// for every lines
	for(unsigned line=0; line<lineCount; line++)
	{
		FileInfo info;
		
		// Files
		info.readingFileName = localToGlobalFileName(dataPath, readingFileNames[line]);
		if (referenceFileNames)
			info.referenceFileName = localToGlobalFileName(dataPath, (*referenceFileNames)[line]);
		if (configFileNames)
			info.configFileName = localToGlobalFileName(configPath, (*configFileNames)[line]);
		
		// Load transformations
		if(found3dInitialTrans)
			info.initialTransformation = getTransform(data, "iT", 3, line);
		if(found2dInitialTrans)
			info.initialTransformation = getTransform(data, "iT", 2, line);
		if(found3dGroundTruthTrans)
			info.groundTruthTransformation = getTransform(data, "gT", 3, line);
		if(found2dGroundTruthTrans)
			info.groundTruthTransformation = getTransform(data, "gT", 2, line);
		
		// Build the list
		this->push_back(info);
	}
	
	// Debug: Print the list
	/*for(unsigned i=0; i<list.size(); i++)
	{
		cout << "\n--------------------------" << endl;
		cout << "Sequence " << i << ":" << endl;
		cout << "Reading path: " << list[i].readingFileName << endl;
		cout << "Reference path: " << list[i].referenceFileName << endl;
		cout << "Extension: " << list[i].fileExtension << endl;
		cout << "Tranformation:\n" << list[i].initialTransformation << endl;
		cout << "Grativity:\n" << list[i].gravity << endl;
	}
	*/
}

//! Join parentPath and fileName and return the result as a global path
template<typename T>
std::string PointMatcherIO<T>::FileInfoVector::localToGlobalFileName(const std::string& parentPath, const std::string& fileName)
{
	std::string globalFileName(fileName);
	if (!boost::filesystem::exists(globalFileName))
	{
		const boost::filesystem::path globalFilePath(boost::filesystem::path(parentPath) /  boost::filesystem::path(fileName));
		#if BOOST_FILESYSTEM_VERSION >= 3
		globalFileName = globalFilePath.string();
		#else
		globalFileName = globalFilePath.file_string();
		#endif
	}
	validateFile(globalFileName);
	return globalFileName;
}

//! Return whether there is a valid transformation named prefix in data
template<typename T>
bool PointMatcherIO<T>::FileInfoVector::findTransform(const CsvElements& data, const std::string& prefix, unsigned dim)
{
	bool found(true);
	for(unsigned i=0; i<dim+1; i++)
	{
		for(unsigned j=0; j<dim+1; j++)
		{
			stringstream transName;
			transName << prefix << i << j;
			found = found && (data.find(transName.str()) != data.end());
		}
	}
	return found;
}

//! Return the transformation named prefix from data
template<typename T>
typename PointMatcherIO<T>::TransformationParameters PointMatcherIO<T>::FileInfoVector::getTransform(const CsvElements& data, const std::string& prefix, unsigned dim, unsigned line)
{
	TransformationParameters transformation(TransformationParameters::Identity(dim+1, dim+1));
	for(unsigned i=0; i<dim+1; i++)
	{
		for(unsigned j=0; j<dim+1; j++)
		{
			stringstream transName;
			transName << prefix << i << j;
			CsvElements::const_iterator colIt(data.find(transName.str()));
			const T value = boost::lexical_cast<T> (colIt->second[line]);
			transformation(i,j) = value;
		}
	}
	return transformation;
}

template struct PointMatcherIO<float>::FileInfoVector;
template struct PointMatcherIO<double>::FileInfoVector;

//! Throw a runtime_error exception if fileName cannot be opened
void PointMatcherSupport::validateFile(const std::string& fileName)
{
	boost::filesystem::path fullPath(fileName);

	ifstream ifs(fileName.c_str());
	if (!ifs.good() || !boost::filesystem::is_regular_file(fullPath))
	#if BOOST_FILESYSTEM_VERSION >= 3
		#if BOOST_VERSION >= 105000
				throw runtime_error(string("Cannot open file ") + boost::filesystem::complete(fullPath).generic_string());
		#else
				throw runtime_error(string("Cannot open file ") + boost::filesystem3::complete(fullPath).generic_string());
		#endif
	#else
		throw runtime_error(string("Cannot open file ") + boost::filesystem::complete(fullPath).native_file_string());
    #endif
}


//! Load a point cloud from a file, determine format from extension
template<typename T>
typename PointMatcher<T>::DataPoints PointMatcher<T>::DataPoints::load(const std::string& fileName)
{
	const boost::filesystem::path path(fileName);
	const string& ext(boost::filesystem::extension(path));
	if (boost::iequals(ext, ".vtk"))
		return PointMatcherIO<T>::loadVTK(fileName);
	else if (boost::iequals(ext, ".csv"))
		return PointMatcherIO<T>::loadCSV(fileName);
	else if (boost::iequals(ext, ".ply"))
		return PointMatcherIO<T>::loadPLY(fileName);
	else if (boost::iequals(ext, ".pcd"))
		return PointMatcherIO<T>::loadPCD(fileName);
	else
		throw runtime_error("loadAnyFormat(): Unknown extension \"" + ext + "\" for file \"" + fileName + "\", extension must be either \".vtk\" or \".csv\"");
}

template
PointMatcher<float>::DataPoints PointMatcher<float>::DataPoints::load(const std::string& fileName);
template
PointMatcher<double>::DataPoints PointMatcher<double>::DataPoints::load(const std::string& fileName);


//! @brief Load comma separated values (csv) file
//! @param fileName a string containing the path and the file name
//! 
//! This loader has 3 behaviors since there is no official standard for
//! csv files. A 2D or 3D point cloud will be created automatically if:
//!   - there is a header with columns named x, y and optionnaly z
//!   - there are only 2 or 3 columns in the file
//!
//! Otherwise, the user is asked to enter column id manually which might 
//! block automatic processing.
template<typename T>
typename PointMatcher<T>::DataPoints PointMatcherIO<T>::loadCSV(const std::string& fileName)
{
	ifstream ifs(fileName.c_str());
	
	validateFile(fileName);

	return loadCSV(ifs);
}

template<typename T>
PointMatcherIO<T>::SupportedLabel::SupportedLabel(const std::string& internalName, const std::string& externalName, const PMPropTypes& type):
	internalName(internalName),
	externalName(externalName),
	type(type)
{
}

template<typename T>
typename PointMatcherIO<T>::SublabelAssociationMap PointMatcherIO<T>::getFeatAssocationMap()
{
	// FIXME: this should be depreciated and replace by SupportedLabel
	const SublabelAssociationMap assoc_map = boost::assign::map_list_of
			("x", LabelAssociationPair(0,"x"))
			("y", LabelAssociationPair(1,"y"))
			("z", LabelAssociationPair(2,"z"))
			("pad", LabelAssociationPair(3,"pad"));
	return assoc_map;
}

template<typename T>
typename PointMatcherIO<T>::SublabelAssociationMap PointMatcherIO<T>::getDescAssocationMap()
{
	// FIXME: this should be depreciated and replace by SupportedLabel
	const SublabelAssociationMap assoc_map = boost::assign::map_list_of
			("nx", LabelAssociationPair(0,"normals"))
			("ny", LabelAssociationPair(1,"normals"))
			("nz", LabelAssociationPair(2,"normals"))
			("normal_x", LabelAssociationPair(0,"normals"))
			("normal_y", LabelAssociationPair(1,"normals"))
			("normal_z", LabelAssociationPair(2,"normals"))
			("densities", LabelAssociationPair(0,"densities"))
			("intensity", LabelAssociationPair(0,"intensity"))
			("red", LabelAssociationPair(0,"color"))
			("green", LabelAssociationPair(1,"color"))
			("blue", LabelAssociationPair(2,"color"))
			("alpha", LabelAssociationPair(3,"color"))
			("eigValues0", LabelAssociationPair(0,"eigValues"))
			("eigValues1", LabelAssociationPair(1,"eigValues"))
			("eigValues2", LabelAssociationPair(2,"eigValues"))
			("eigVectors0X", LabelAssociationPair(0,"eigVectors"))
			("eigVectors0Y", LabelAssociationPair(1,"eigVectors"))
			("eigVectors0Z",LabelAssociationPair(2,"eigVectors"))
			("eigVectors1X", LabelAssociationPair(3,"eigVectors"))
			("eigVectors1Y", LabelAssociationPair(4,"eigVectors"))
			("eigVectors1Z",LabelAssociationPair(5,"eigVectors"))
			("eigVectors2X", LabelAssociationPair(6,"eigVectors"))
			("eigVectors2Y", LabelAssociationPair(7,"eigVectors"))
			("eigVectors2Z",LabelAssociationPair(8,"eigVectors"))
			("normals", LabelAssociationPair(0,"normals"))
			("eigValues", LabelAssociationPair(0,"eigValues"))
			("eigVectors", LabelAssociationPair(0,"eigVectors"))
			("color", LabelAssociationPair(0,"color"));
	return assoc_map;
}

template <typename T>
bool PointMatcherIO<T>::featSublabelRegistered(const std::string& externalName)
{
	return getFeatAssocationMap().count(externalName) > 0;
}

template <typename T>
bool PointMatcherIO<T>::descSublabelRegistered(const std::string& externalName)
{
	return getDescAssocationMap().count(externalName) > 0;
}

template <typename T>
typename PointMatcherIO<T>::LabelAssociationPair PointMatcherIO<T>::getFeatAssociationPair(const std::string& externalName)
{
	return getFeatAssocationMap().find(externalName)->second;
}

template <typename T>
typename PointMatcherIO<T>::LabelAssociationPair PointMatcherIO<T>::getDescAssociationPair(const std::string& externalName)
{
	return getDescAssocationMap().find(externalName)->second;
}

template<typename T>
typename PointMatcherIO<T>::PMPropTypes PointMatcherIO<T>::getPMType(const std::string& externalName)
{
	if (featSublabelRegistered(externalName))
		return FEATURE;
	else if (descSublabelRegistered(externalName))
		return DESCRIPTOR;
	else
		return UNSUPPORTED;
	//TODO: add time here
}

// Class LabelGenerator
template<typename T>
void PointMatcherIO<T>::LabelGenerator::add(std::string internalName)
{
	bool findLabel = false;
	for(size_t i=0; i<labels.size(); ++i)
	{
		if(internalName == labels[i].text)
		{
			labels[i].span++;
			findLabel = true;
			break;
		}
		
	}

	if(!findLabel)
	{
		labels.push_back(Label(internalName,1));
	}
}

// Class LabelGenerator
template<typename T>
typename PointMatcher<T>::DataPoints::Labels PointMatcherIO<T>::LabelGenerator::getLabels() const
{
	return labels;
}

template
class PointMatcherIO<float>::LabelGenerator;
template
class PointMatcherIO<double>::LabelGenerator;
template <typename T>


std::string PointMatcherIO<T>::getColLabel(const Label& label, const int row)
{
	std::string externalName;
	if (label.text == "normals")
	{
		if (row == 0)
		{
			externalName = "nx";
		}
		if (row == 1)
		{
			externalName = "ny";
		}
		if (row == 2)
		{
			externalName = "nz";
		}
	}
	else if (label.text == "color")
	{
		if (row == 0)
		{
			externalName = "red";
		}
		if (row == 1)
		{
			externalName = "green";
		}
		if (row == 2)
		{
			externalName = "blue";
		}
		if (row == 3)
			externalName = "alpha";
	}
	else if (label.text == "eigValues")
	{
		externalName = "eigValues" + boost::lexical_cast<string>(row);
	}
	else if (label.text == "eigVectors")
	{
		// format: eigVectors<0-2><X-Z>
		externalName = "eigVectors" + boost::lexical_cast<string>(row/3);

		int row_mod = row % 3;
		if (row_mod == 0)
			externalName += "X";
		else if (row_mod == 1)
			externalName += "Y";
		else if (row_mod == 2)
			externalName += "Z";
	}
	else if (label.span  == 1)
	{
		externalName = label.text;
	}
	else
		externalName = label.text + boost::lexical_cast<std::string>(row);

	return externalName;
}


//! @brief Load comma separated values (csv) file
//! @see loadCSV()
template<typename T>
typename PointMatcher<T>::DataPoints PointMatcherIO<T>::loadCSV(std::istream& is)
{
	vector<GenericInputHeader> csvHeader;
	LabelGenerator featLabelGen, descLabelGen, timeLabelGen;
	Matrix features;
	Matrix descriptors;
	Int64Matrix times;
	unsigned int csvCol = 0;
	unsigned int csvRow = 0;
	bool hasHeader(false);
	bool firstLine(true);
	

	//count lines in the file
	is.unsetf(std::ios_base::skipws);
	unsigned int line_count = std::count(
	        std::istream_iterator<char>(is),
			std::istream_iterator<char>(), 
			'\n');
	//reset the stream
	is.clear();
	is.seekg(0, ios::beg);

	char delimiters[] = " \t,;";
	char *token;
	while (!is.eof())
	{
		string line;
		safeGetLine(is, line);
		// Skip empty lines
		if(line.empty())
			break;
	
		// Look for text header
		unsigned int len = strspn(line.c_str(), " ,+-.1234567890Ee");
		if(len != line.length())
		{
			//cout << "Header detected" << endl;
			hasHeader = true;
		}
		else
		{
			hasHeader = false;
		}

		// Count dimension using first line
		if(firstLine)
		{
			unsigned int dim = 0;
			char tmpLine[1024]; //FIXME: might be problematic for large file
			strcpy(tmpLine, line.c_str());
			char *brkt = 0;
			token = strtok_r(tmpLine, delimiters, &brkt);

			//1- BUILD HEADER
			while (token)
			{
				// Load text header
				if(hasHeader)
				{
					csvHeader.push_back(GenericInputHeader(string(token)));
				}
				dim++;
				token = strtok_r(NULL, delimiters, &brkt);
			}
			
			if (!hasHeader)
			{
				// Check if it is a simple file with only coordinates
				if (!(dim == 2 || dim == 3))
				{
					int idX=0, idY=0, idZ=0;

					cout << "WARNING: " << dim << " columns detected. Not obvious which columns to load for x, y or z." << endl;
					cout << endl << "Enter column ID (starting from 0) for x: ";
					cin >> idX;
					cout << "Enter column ID (starting from 0) for y: ";
					cin >> idY;
					cout << "Enter column ID (starting from 0, -1 if 2D data) for z: ";
					cin >> idZ;

					// Fill with unkown column names
					for(unsigned int i=0; i<dim; i++)
					{
						std::ostringstream os;
						os << "empty" << i;

						csvHeader.push_back(GenericInputHeader(os.str()));
					}
					
					// Overwrite with user inputs
					csvHeader[idX] = GenericInputHeader("x");
					csvHeader[idY] = GenericInputHeader("y");
					if(idZ != -1)
						csvHeader[idZ] = GenericInputHeader("z");
				}
				else
				{
					// Assume logical order...
					csvHeader.push_back(GenericInputHeader("x"));
					csvHeader.push_back(GenericInputHeader("y"));
					if(dim == 3)
						csvHeader.push_back(GenericInputHeader("z"));
				}
			}

			//2- PROCESS HEADER
			// Load known features, descriptors, and time
			const SupportedLabels externalLabels = getSupportedExternalLabels();

			// Counters
			int rowIdFeatures = 0;
			int rowIdDescriptors = 0;
			int rowIdTime = 0;

			
			// Loop through all known external names (ordered list)
			for(size_t i=0; i<externalLabels.size(); i++)
			{
				const SupportedLabel supLabel = externalLabels[i];

				for(size_t j=0; j < csvHeader.size(); j++)
				{
					if(supLabel.externalName == csvHeader[j].name)
					{
						csvHeader[j].isKnownName = true;
						csvHeader[j].matrixType = supLabel.type;

						switch (supLabel.type)
						{
							case FEATURE:
								csvHeader[j].matrixRowId = rowIdFeatures;
								featLabelGen.add(supLabel.internalName);
								rowIdFeatures++;
								break;
							case DESCRIPTOR:
								csvHeader[j].matrixRowId = rowIdDescriptors;
								descLabelGen.add(supLabel.internalName);
								rowIdDescriptors++;
								break;
							case TIME:
								csvHeader[j].matrixRowId = rowIdTime;
								timeLabelGen.add(supLabel.internalName);
								rowIdTime++;
								break;
							default:
								throw runtime_error(string("CSV parse error: encounter a type different from FEATURE, DESCRIPTOR and TIME. Implementation not supported. See the definition of 'enum PMPropTypes'"));
								break;
						}
					}
				}
			}

			// loop through the remaining UNSUPPORTED labels and assigned them to a descriptor row
			for(unsigned int i=0; i<csvHeader.size(); i++)
			{
				if(csvHeader[i].matrixType == UNSUPPORTED)
				{
					csvHeader[i].matrixRowId = rowIdDescriptors;
					csvHeader[i].matrixType = DESCRIPTOR;
					descLabelGen.add(csvHeader[i].name);
					rowIdDescriptors++;
				}
			}
		

			//3- RESERVE MEMORY
			if(hasHeader && line_count > 0)
				line_count--;

			const unsigned int featDim = featLabelGen.getLabels().totalDim();
			const unsigned int descDim = descLabelGen.getLabels().totalDim();
			const unsigned int timeDim = timeLabelGen.getLabels().totalDim();
			const unsigned int nbPoints = line_count;

			features = Matrix(featDim, nbPoints);
			descriptors = Matrix(descDim, nbPoints);
			times = Int64Matrix(timeDim, nbPoints);
		}


		//4- LOAD DATA (this start again from the first line)
		char* brkt = 0;
		char line_c[1024];//FIXME: this might be a problem for large files
		strcpy(line_c,line.c_str());
		token = strtok_r(line_c, delimiters, &brkt);
		
		if(!(hasHeader && firstLine))
		{
			// Parse a line
			csvCol = 0;
			while (token)
			{
				if(csvCol > (csvHeader.size() - 1))
				{
					// Error check (too much data)
					throw runtime_error(
					(boost::format("CSV parse error: at line %1%, too many elements to parse compare to the header number of columns (col=%2%).") % csvRow % csvHeader.size()).str());
				}
				
				// Alias
				const int matrixRow = csvHeader[csvCol].matrixRowId;
				const int matrixCol = csvRow;
				
				switch (csvHeader[csvCol].matrixType)
				{
					case FEATURE:
						features(matrixRow, matrixCol) = lexical_cast_scalar_to_string<T>(string(token));
						break;
					case DESCRIPTOR:
						descriptors(matrixRow, matrixCol) = lexical_cast_scalar_to_string<T>(token);
						break;
					case TIME:
						times(matrixRow, matrixCol) = lexical_cast_scalar_to_string<boost::int64_t>(token);
						break;
					default:
						throw runtime_error(string("CSV parse error: encounter a type different from FEATURE, DESCRIPTOR and TIME. Implementation not supported. See the definition of 'enum PMPropTypes'"));
						break;

				}

				//fetch next element
				token = strtok_r(NULL, delimiters, &brkt);
				csvCol++;
			}
		
		
			// Error check (not enough data)
			if(csvCol != (csvHeader.size()))
			{
				throw runtime_error(
				(boost::format("CSV parse error: at line %1%, not enough elements to parse compare to the header number of columns (col=%2%).") % csvRow % csvHeader.size()).str());
			}

			csvRow++;
		}

		firstLine = false;
		
	}

	// 5- ASSEMBLE FINAL DATAPOINTS
	DataPoints loadedPoints(features, featLabelGen.getLabels());

	if (descriptors.rows() > 0)
	{
		loadedPoints.descriptors = descriptors;
		loadedPoints.descriptorLabels = descLabelGen.getLabels();
	}

	if(times.rows() > 0)
	{
		loadedPoints.times = times;
		loadedPoints.timeLabels = timeLabelGen.getLabels();	
	}

	// Ensure homogeous coordinates
	if(!loadedPoints.featureExists("pad"))
	{
		loadedPoints.addFeature("pad", Matrix::Ones(1,features.cols()));
	}

	return loadedPoints;
}

template
PointMatcher<float>::DataPoints PointMatcherIO<float>::loadCSV(const std::string& fileName);
template
PointMatcher<double>::DataPoints PointMatcherIO<double>::loadCSV(const std::string& fileName);

//! Save a point cloud to a file, determine format from extension
template<typename T>
void PointMatcher<T>::DataPoints::save(const std::string& fileName, bool binary) const
{
	const boost::filesystem::path path(fileName);
	const string& ext(boost::filesystem::extension(path));
	if (boost::iequals(ext, ".vtk"))
		return PointMatcherIO<T>::saveVTK(*this, fileName, binary);

	if (binary)
		throw runtime_error("save(): Binary writing is not supported together with extension \"" + ext + "\". Currently binary writing is only supported with \".vtk\".");

	if (boost::iequals(ext, ".csv"))
		return PointMatcherIO<T>::saveCSV(*this, fileName);
	else if (boost::iequals(ext, ".ply"))
		return PointMatcherIO<T>::savePLY(*this, fileName);
	else if (boost::iequals(ext, ".pcd"))
		return PointMatcherIO<T>::savePCD(*this, fileName);
	else
		throw runtime_error("save(): Unknown extension \"" + ext + "\" for file \"" + fileName + "\", extension must be either \".vtk\", \".ply\", \".pcd\" or \".csv\"");
}

template
void PointMatcher<float>::DataPoints::save(const std::string& fileName, bool binary) const;
template
void PointMatcher<double>::DataPoints::save(const std::string& fileName, bool binary) const;

//! Save point cloud to a file as CSV
template<typename T>
void PointMatcherIO<T>::saveCSV(const DataPoints& data, const std::string& fileName)
{
	ofstream ofs(fileName.c_str());
	if (!ofs.good())
		throw runtime_error(string("Cannot open file ") + fileName);
	saveCSV(data, ofs);
}

//! Save point cloud to a stream as CSV
template<typename T>
void PointMatcherIO<T>::saveCSV(const DataPoints& data, std::ostream& os)
{
	const int pointCount(data.features.cols());
	const int dimCount(data.features.rows());
	const int descDimCount(data.descriptors.rows());
	
	if (pointCount == 0)
	{
		cerr << "Warning, no points, doing nothing" << endl;
		return;
	}
	
	// write header
	for (int i = 0; i < dimCount - 1; i++)
	{
		os << data.featureLabels[i].text;

		if (!((i == (dimCount - 2)) && descDimCount == 0))
			os << ",";
	}

	int n = 0;
	for (size_t i = 0; i < data.descriptorLabels.size(); i++)
	{
		Label lab = data.descriptorLabels[i];
		for (size_t s = 0; s < lab.span; s++)
		{
			os << getColLabel(lab,s);
			if (n != (descDimCount - 1))
				os << ",";
			n++;
		}
	}

	os << "\n";

	// write points
	for (int p = 0; p < pointCount; ++p)
	{
		for (int i = 0; i < dimCount-1; ++i)
		{
			os << data.features(i, p);
			if(!((i == (dimCount - 2)) && descDimCount == 0))
				os << " , ";
		}

		for (int i = 0; i < descDimCount; i++)
		{
			os << data.descriptors(i,p);
			if (i != (descDimCount - 1))
				os << ",";
		}
		os << "\n";
	}
	
}

template
void PointMatcherIO<float>::saveCSV(const DataPoints& data, const std::string& fileName);
template
void PointMatcherIO<double>::saveCSV(const DataPoints& data, const std::string& fileName);

//! Load point cloud from a file as VTK
template<typename T>
typename PointMatcher<T>::DataPoints PointMatcherIO<T>::loadVTK(const std::string& fileName)
{
	ifstream ifs(fileName.c_str());
	if (!ifs.good())
		throw runtime_error(string("Cannot open file ") + fileName);
	return loadVTK(ifs);
}

void skipBlock(bool binary, int binarySize, std::istream & is, bool hasSeparateSizeParameter = true){
	int n;
	int size;
	is >> n;
	if(hasSeparateSizeParameter) {
		is >> size;
	} else {
		size = n;
	}

	std::string line;
	getline(is, line); // remove line end after parameters;
	if(binary){
		is.seekg(size * binarySize, std::ios_base::cur);
	} else {
		for (int p = 0; p < n; p++)
		{
			getline(is, line);
		}
	}
}

//! Load point cloud from a stream as VTK
template<typename T>
typename PointMatcher<T>::DataPoints PointMatcherIO<T>::loadVTK(std::istream& is)
{
	//typedef typename DataPoints::Label Label;
	//typedef typename DataPoints::Labels Labels;
	
	DataPoints loadedPoints;

	// parse header
	string line;
	getline(is, line);
	if (line.find("# vtk DataFile Version") != 0)
		throw runtime_error(string("Wrong magic header, found ") + line);
	getline(is, line);
	getline(is, line);

	const bool isBinary = (line == "BINARY");
	if (line != "ASCII"){
		if(!isBinary){
			throw runtime_error(string("Wrong file type, expecting ASCII or BINARY, found ") + line);
		}
	}
	getline(is, line);

	SupportedVTKDataTypes dataType;
	if (line == "DATASET POLYDATA")
		dataType = POLYDATA;
	else if (line == "DATASET UNSTRUCTURED_GRID")
		dataType = UNSTRUCTURED_GRID;
	else
		throw runtime_error(string("Wrong data type, expecting DATASET POLYDATA, found ") + line);


	// parse points and descriptors
	string fieldName;
	string name;
	int dim = 0;
	int pointCount = 0;
	string type;
	while (is.good())
	{
		is >> fieldName;
		
		// load features
		if(fieldName == "POINTS")
		{
			is >> pointCount;
			is >> type;
			getline(is, line); // remove line end after parameters!

			if(!(type == "float" || type == "double"))
					throw runtime_error(string("Field POINTS can only be of type double or float"));

			Matrix features(4, pointCount);
			for (int p = 0; p < pointCount; ++p)
			{
				readVtkData(type, isBinary, features.template block<3, 1>(0, p), is);
				features(3, p) = 1.0;
			}
			loadedPoints.addFeature("x", features.row(0));
			loadedPoints.addFeature("y", features.row(1));
			loadedPoints.addFeature("z", features.row(2));
			loadedPoints.addFeature("pad", features.row(3));
		}

		//////////////////////////////////////////////////////////
		// Dataset type
		// POLYDATA
		else if(dataType == POLYDATA && fieldName == "VERTICES")
		{
			skipBlock(isBinary, 4, is);
		}

		else if(dataType == POLYDATA && fieldName == "LINES")
		{
			skipBlock(isBinary, 4, is);
		}

		else if(dataType == POLYDATA && fieldName == "POLYGONS")
		{
			skipBlock(isBinary, 4, is);
		}

		else if(dataType == POLYDATA && fieldName == "TRIANGLE_STRIPS")
		{
			skipBlock(isBinary, 4, is);
		}

		// Unstructure Grid
		else if(dataType == UNSTRUCTURED_GRID && fieldName == "CELLS")
		{
			skipBlock(isBinary, 4, is);
		}
		else if(dataType == UNSTRUCTURED_GRID && fieldName == "CELL_TYPES")
		{
			skipBlock(isBinary, 4, is, false); // according to http://www.vtk.org/VTK/img/file-formats.pdf CELL_TYPES only has one parameter (n)
		}

		//////////////////////////////////////////////////////////
		// Point data
		else if(fieldName == "POINT_DATA")
		{
			int descriptorCount;
			is >> descriptorCount;
			if(pointCount != descriptorCount)
				throw runtime_error(string("The size of POINTS is different than POINT_DATA"));
		}
		//////////////////////////////////////////////////////////
		// Field data is ignored
		else if (fieldName == "FIELD")
		{
			string fieldDataName;
			int fieldDataCount;
			is >> fieldDataName >> fieldDataCount;

			for (int f = 0; f < fieldDataCount; f++)
			{
				//getline(is, line);
				int numTuples;
				is >> name >> dim >> numTuples >> type;

				if(type == "vtkIdType") // skip that type
				{
					if(isBinary){
						is.seekg(dim * numTuples * 4, std::ios_base::cur);
					} else {
						int t_val;
						for (int t = 0; t < dim * numTuples; t++ )
						{
							is >> t_val;
						}
					}
				}
				else if(!(type == "float" || type == "double"))
						throw runtime_error(string("Field " + fieldName + " is " + type + " but can only be of type double or float"));
						 

				Matrix descriptor(dim, pointCount);
				readVtkData(type, isBinary, descriptor.transpose(), is);
				loadedPoints.addDescriptor(name, descriptor);
			}
		}
		else // Load descriptors
		{
			// descriptor name
			is >> name;

			bool skipLookupTable = false;
			bool isColorScalars = false;
			if(fieldName == "SCALARS")
			{
				dim = 1;
				is >> type;
				skipLookupTable = true;
			}
			else if(fieldName == "VECTORS")
			{
				dim = 3;
				is >> type;
			}
			else if(fieldName == "TENSORS")
			{
				dim = 9;
				is >> type;
			}
			else if(fieldName == "NORMALS")
			{
				dim = 3;
				is >> type;
			}
			else if(fieldName == "COLOR_SCALARS")
			{
				is >> dim;
				type = "float";
				isColorScalars = true;
			}
			else
				throw runtime_error(string("Unknown field name " + fieldName + ", expecting SCALARS, VECTORS, TENSORS, NORMALS or COLOR_SCALARS."));

			
			getline(is, line); // remove rest of the parameter line including its line end;

			Matrix descriptor(dim, pointCount);
			if(isColorScalars && isBinary) {
				std::vector<unsigned char> buffer(dim);
				for (int i = 0; i < pointCount; ++i){
					is.read(reinterpret_cast<char *>(&buffer.front()), dim);
					for(int r=0; r < dim; ++r){
						descriptor(r, i) = buffer[r] / static_cast<T>(255.0);
					}
				}
			} else {
				if(!(type == "float" || type == "double"))
						throw runtime_error(string("Field " + fieldName + " is " + type + " but can only be of type double or float"));

				// Skip LOOKUP_TABLE line
				if(skipLookupTable)
				{
					getline(is, line);
				}
				readVtkData(type, isBinary, descriptor.transpose(), is);
			}
			loadedPoints.addDescriptor(name, descriptor);
		}
	}
	
	return loadedPoints;
}

template
PointMatcherIO<float>::DataPoints PointMatcherIO<float>::loadVTK(const std::string& fileName);
template
PointMatcherIO<double>::DataPoints PointMatcherIO<double>::loadVTK(const std::string& fileName);


//! Save point cloud to a file as VTK
template<typename T>
void PointMatcherIO<T>::saveVTK(const DataPoints& data, const std::string& fileName, bool binary)
{
	typedef typename InspectorsImpl<T>::VTKFileInspector VTKInspector;
	
	Parametrizable::Parameters param;
	boost::assign::insert(param) ("baseFileName", "");
	boost::assign::insert(param) ("writeBinary", toParam(binary));
	VTKInspector vtkInspector(param);
	vtkInspector.dumpDataPoints(data, fileName);
}


template
void PointMatcherIO<float>::saveVTK(const PointMatcherIO<float>::DataPoints& data, const std::string& fileName, bool binary);
template
void PointMatcherIO<double>::saveVTK(const PointMatcher<double>::DataPoints& data, const std::string& fileName, bool binary);

//! @brief Load polygon file format (ply) file
//! @param fileName a string containing the path and the file name
//!
//! Note: that the PLY does not define a standard for point clouds
//! Warning: Binary PLY files are not supported, only ASCII
//! Only PLY files with elements named "vertex" are supported
//! "vertex" should have 2 or 3 properties names "x", "y", "z" to define features.
//!
template<typename T>
typename PointMatcherIO<T>::DataPoints PointMatcherIO<T>::loadPLY(const std::string& fileName)
{
	ifstream ifs(fileName.c_str());
	if (!ifs.good())
		throw runtime_error(string("Cannot open file ") + fileName);
	return loadPLY(ifs);
}

template
PointMatcherIO<float>::DataPoints PointMatcherIO<float>::loadPLY(const string& fileName);
template
PointMatcherIO<double>::DataPoints PointMatcherIO<double>::loadPLY(const string& fileName);

//! @brief Load polygon file format (PLY) file
//! @see loadPLY()
template <typename T>
typename PointMatcherIO<T>::DataPoints PointMatcherIO<T>::loadPLY(std::istream& is)
{
	//TODO: adapt following loadCSV()
	typedef vector<PLYElement*> Elements;

	/*
	Steps:
	1- PARSE PLY HEADER
	2- ASSIGN PLY PROPERTIES TO DATAPOINTS ROWS
	3- Reserve memory for a DataPoints
	4- Parse PLY vertex to appropriate DataPoints cols and rows 
	5- Assemble final DataPoints

	PLY organisation:

	element 1 [name, size]
		- property 1 [type, name]
		- property 2
		- ...
	element 2
		-...

	*/


	///////////////////////////
	// 1- PARSE PLY HEADER
	bool format_defined = false;
	bool header_processed = false;

	Elements elements;
	PLYElementF element_f; // factory
	PLYElement* current_element = NULL;
	bool skip_props = false; // flag to skip properties if element is not supported
	unsigned elem_offset = 0; // keep track of line position of elements that are supported
	string line;
	getline(is, line);

	if (line.find("ply") != 0) {
		throw runtime_error(string("PLY parse error: wrong magic header, found <") + line + string(">"));
	}

	while (!header_processed)
	{
		if(!getline(is, line))
			throw runtime_error("PLY parse error: reached end of file before end of header definition");


		if ( line.empty() )
			continue;
		istringstream stringstream (line);

		string keyword;
		stringstream >> keyword;

		// ignore comment
		if (keyword == "comment") {
			continue;
		}

		// We only support ascii and ply version 1.0
		if (keyword == "format")
		{
			if (format_defined)
				throw runtime_error("PLY parse error: format already defined");

			string format_str, version_str;
			stringstream >> format_str >> version_str;

			if (format_str != "ascii" && format_str != "binary_little_endian" && format_str != "binary_big_endian")
				throw runtime_error(string("PLY parse error: format <") + format_str + string("> is not supported"));

			if (format_str == "binary_little_endian" || format_str == "binary_big_endian")
				throw runtime_error(string("PLY parse error: binary PLY files are not supported"));
			if (version_str != "1.0")
			{
				throw runtime_error(string("PLY parse error: version <") + version_str + string("> of ply is not supported"));
			}

			format_defined = true;

		}
		else if (keyword == "element")
		{
			

			string elem_name, elem_num_s;
			stringstream >> elem_name >> elem_num_s;

			unsigned elem_num;
			try
			{
				elem_num = boost::lexical_cast<unsigned>(elem_num_s);
			}
			catch (boost::bad_lexical_cast& e)
			{
				throw runtime_error(string("PLY parse error: bad number of elements ") + elem_num_s + string(" for element ") + elem_name);
			}

			if (element_f.elementSupported(elem_name))
			{
				// Initialize current element
				PLYElement* elem = element_f.createElement(elem_name, elem_num, elem_offset);
				current_element = elem;

				// check that element is not already defined
				for (typename Elements::const_iterator it = elements.begin(); it != elements.end(); it++ )
				{
					if (**it == *elem) {
						throw runtime_error(string("PLY parse error: element: ") + elem_name + string( "is already defined"));
					}
				}
				elements.push_back(elem);
				skip_props = false;
			}
			else
			{
				LOG_WARNING_STREAM("PLY parse warning: element " << elem_name << " not supported. Skipping.");
				skip_props = true;
			}

			elem_offset += elem_num;
		}
		else if (keyword == "property")
		{
			if (current_element == NULL)
			{
				throw runtime_error("PLY parse error: property listed without defining an element");
			}

			if (skip_props)
				continue;

			string next, prop_type, prop_name;
			stringstream >> next;

			// PLY list property
			if (next == "list")
			{
				string prop_idx_type;
				stringstream >> prop_idx_type >> prop_type >> prop_name;

				PLYProperty list_prop(prop_idx_type, prop_type, prop_name, current_element->total_props);

				current_element->addProperty(list_prop);
			}
			// PLY regular property
			else
			{
				prop_type = next;
				stringstream >> prop_name;
				PLYProperty prop(prop_type, prop_name, current_element->total_props);

				current_element->addProperty(prop);
			}

			current_element->total_props++;
		}
		else if (keyword == "end_header")
		{
			if (!format_defined)
			{
				throw runtime_error(string("PLY parse error: format not defined in header"));
			}

			if (elements.size() == 0)
			{
				throw runtime_error(string("PLY parse error: no elements defined in header"));
			}

			header_processed = true;
		}
	}

	///////////////////////////
	// 2- ASSIGN PLY PROPERTIES TO DATAPOINTS ROWS
	
	// Fetch vertex
	PLYElement* vertex = elements[0];
	
	if(vertex->name != "vertex")
	{
		throw runtime_error(string("PLY parse error: vertex should be the first element defined."));
	}
		
	// Known features and descriptors
	const SupportedLabels & externalLabels = getSupportedExternalLabels();
	
	int rowIdFeatures = 0;
	int rowIdDescriptors = 0;
	
	LabelGenerator featLabelGen, descLabelGen;
	
	// Loop through all known external names (ordered list)
	for(size_t i=0; i<externalLabels.size(); i++)
	{
		const SupportedLabel & supLabel = externalLabels[i];

		//Search if that feature exist
		for(it_PLYProp it=vertex->properties.begin(); it!=vertex->properties.end(); ++it)
		{
			if(supLabel.externalName == it->name)
			{
				// Assign rowId in that order
				if(supLabel.type == FEATURE)
				{
					it->pmRowID = rowIdFeatures;

					// Prepare feature labels
					featLabelGen.add(supLabel.internalName);

					rowIdFeatures++;
				}
				else if (supLabel.type == DESCRIPTOR)
				{
					it->pmRowID = rowIdDescriptors;

					// Prepare descriptor labels
					descLabelGen.add(supLabel.internalName);

					rowIdDescriptors++;
				}
				else
				{
					throw runtime_error(string("PLY parse error: encounter a type different from FEATURE and DESCRIPTOR. Implementation not supported. See the definition of 'enum PMPropTypes'"));
				}
				break;
			}
		}

		//TODO: Handle random descriptor names
	}

	///////////////////////////
	// 3- RESERVE DATAPOINTS MEMORY

	const int featDim = featLabelGen.getLabels().totalDim();
	const int descDim = descLabelGen.getLabels().totalDim();
	const int nbPoints = vertex->num;

	Matrix features = Matrix(featDim, nbPoints);
	Matrix descriptors = Matrix(descDim, nbPoints);
	//TODO: add time



	///////////////////////////
	// 4- PARSE PLY DATA (vertex)
	const int nbProp = vertex->total_props;
	const int nbValues = nbPoints*nbProp;
	int propID = 0;
	int col = 0;
	for(int i=0; i<nbValues; i++)
	{
		T value;
		if(!(is >> value))
		{
			throw runtime_error(
			(boost::format("PLY parse error: expected %1% values (%2% points with %3% properties) but only found %4% values.") % nbValues % nbPoints % nbProp % i).str());
		}
		else
		{
			const int row = vertex->properties[propID].pmRowID;
			const PMPropTypes type = vertex->properties[propID].pmType;
			
			if (vertex->properties[propID].name == "red" || vertex->properties[propID].name == "green" || vertex->properties[propID].name == "blue" || vertex->properties[propID].name == "alpha") {
				value /= 255.0;
			}

			if(type == FEATURE)
			{
				features(row, col) = value;
			}
			else if(type == DESCRIPTOR)
			{
				descriptors(row, col) = value;
			}

			++propID;

			if(propID >= nbProp)
			{
				propID = 0;
				++col;
			}
		}
	}



	///////////////////////////
	// 5- ASSEMBLE FINAL DATAPOINTS
	
	DataPoints loadedPoints;

	if (descriptors.rows() > 0)
	{
		loadedPoints = DataPoints(features, featLabelGen.getLabels(), 
		                          descriptors,descLabelGen.getLabels());
	}
	else
	{
		DataPoints loadedPoints(features, featLabelGen.getLabels());
	}

	// Ensure homogeous coordinates
	if(!loadedPoints.featureExists("pad"))
	{
		loadedPoints.addFeature("pad", Matrix::Ones(1,nbPoints));
	}

	return loadedPoints;

}

template<typename T>
void PointMatcherIO<T>::savePLY(const DataPoints& data,
		const std::string& fileName)
{
	//typedef typename DataPoints::Labels Labels;

	ofstream ofs(fileName.c_str());
	if (!ofs.good())
		throw runtime_error(string("Cannot open file ") + fileName);

	const int pointCount(data.features.cols());
	const int featCount(data.features.rows());
	const int descRows(data.descriptors.rows());


	if (pointCount == 0)
	{
		cerr << "Warning, no points, doing nothing" << endl;
		return;
	}

	ofs << "ply\n" <<"format ascii 1.0\n";
	ofs << "element vertex " << pointCount << "\n";
	for (int f=0; f <(featCount-1); f++)
	{
		ofs << "property float " << data.featureLabels[f].text << "\n";
	}

	for (size_t i = 0; i < data.descriptorLabels.size(); i++)
	{
		Label lab = data.descriptorLabels[i];
		for (size_t s = 0; s < lab.span; s++)
		{
			ofs << "property float " << getColLabel(lab,s) << "\n";
		}
	}

	ofs << "end_header\n";

	// write points
	for (int p = 0; p < pointCount; ++p)
	{
		for (int f = 0; f < featCount - 1; ++f)
		{
			ofs << data.features(f, p);
			if(!(f == featCount-2 && descRows == 0))
				ofs << " ";
		}

		bool datawithColor = data.descriptorExists("color");
		int colorStartingRow = data.getDescriptorStartingRow("color");
		int colorEndRow = colorStartingRow + data.getDescriptorDimension("color");
		for (int d = 0; d < descRows; ++d)
		{
			if (datawithColor && d >= colorStartingRow && d < colorEndRow) {
				ofs << static_cast<unsigned>(data.descriptors(d, p) * 255.0);
			} else {
				ofs << data.descriptors(d, p);
			}
			if(d != descRows-1)
				ofs << " ";
		}
		ofs << "\n";
	}

	ofs.close();
}

template
void PointMatcherIO<float>::savePLY(const DataPoints& data, const std::string& fileName);
template
void PointMatcherIO<double>::savePLY(const DataPoints& data, const std::string& fileName);

//! @(brief) Regular PLY property constructor
template<typename T>
PointMatcherIO<T>::PLYProperty::PLYProperty(const std::string& type,
		const std::string& name, const unsigned pos, const bool is_feature) :
		name(name), 
		type(type), 
		pos(pos), 
		is_feature(is_feature)  
{
	if (plyPropTypeValid(type))
	{
		is_list = false;
	}
	else
	{
		throw std::runtime_error(
				std::string("PLY parse error: property type ") + type
						+ std::string(" for property ") + name
						+ std::string(" is invalid"));
	}

	pmType = getPMType(name);
	pmRowID = -1;

}

//! @(brief) PLY property list constructor
template<typename T>
PointMatcherIO<T>::PLYProperty::PLYProperty(const std::string& idx_type,
		const std::string& type, const std::string& name, const unsigned pos, const bool is_feature) :
		name(name), 
		type(type), 
		idx_type(idx_type), 
		pos(pos), 
		is_feature(is_feature)
{
	if (plyPropTypeValid(idx_type) && plyPropTypeValid(type)) 
	{
		is_list = true;
	} 
	else
	{
		throw std::runtime_error(
				std::string("PLY parse error: property list type ") + idx_type
						+ std::string(" ") + type
						+ std::string(" for property ") + name
						+ std::string(" is invalid"));
	}

	pmType = getPMType(name);
	pmRowID = -1;
}

template
class PointMatcherIO<float>::PLYElement;
template
class PointMatcherIO<double>::PLYElement;

template
class PointMatcherIO<float>::PLYProperty;
template
class PointMatcherIO<double>::PLYProperty;

template <typename T>
void PointMatcherIO<T>::PLYElement::addProperty(
		PLYProperty& prop) 
{
	if (prop.pmType == FEATURE)
	{
		nbFeatures++;
	}
	else if (prop.pmType == DESCRIPTOR)
	{
		nbDescriptors++;
	}
		
	properties.push_back(prop);
}


template <typename T>
bool PointMatcherIO<T>::PLYElement::supportsProperty(const PLYProperty& prop) const
{
	return getPMType(prop.name) != UNSUPPORTED;
}



template <typename T>
typename PointMatcherIO<T>::PLYElementF::ElementTypes PointMatcherIO<T>::PLYElementF::getElementType(const std::string& elem_name)
{
	string lc = elem_name;
	boost::algorithm::to_lower(lc);
	if (lc == "vertex")
	{
		return VERTEX;
	}
	else
	{
		return UNSUPPORTED;
	}
}

template <typename T>
bool PointMatcherIO<T>::PLYElementF::elementSupported(const std::string& elem_name)
{
	return getElementType(elem_name) != UNSUPPORTED;
}

template<typename T>
typename PointMatcherIO<T>::PLYElement* PointMatcherIO<T>::PLYElementF::createElement(
		const std::string& elem_name, const int elem_num, const unsigned offset) {
	ElementTypes type = getElementType(elem_name);
	if (type == VERTEX)
		return new PLYVertex(elem_num, offset);
	else
		return NULL;
}

template<typename T>
bool PointMatcherIO<T>::plyPropTypeValid(const std::string& type) {
	return (type == "char" || type == "uchar" || type == "short"
			|| type == "ushort" || type == "int" || type == "uint"
			|| type == "float" || type == "double");
}

template <typename T>
bool PointMatcherIO<T>::PLYElement::operator==(const PLYElement& rhs) const
{
	return name == rhs.name;
}


template <typename T>
bool PointMatcherIO<T>::PLYProperty::operator==(const PLYProperty& rhs) const
{
	return name == rhs.name && type == rhs.type;
}

//! @brief Load Point Cloud Library (pcd) file
//! @param fileName a string containing the path and the file name
template<typename T>
typename PointMatcherIO<T>::DataPoints PointMatcherIO<T>::loadPCD(const string& fileName) {
	ifstream ifs(fileName.c_str());
	if (!ifs.good())
		throw runtime_error(string("Cannot open file ") + fileName);
	return loadPCD(ifs);
}

template
PointMatcherIO<float>::DataPoints PointMatcherIO<float>::loadPCD(const string& fileName);
template
PointMatcherIO<double>::DataPoints PointMatcherIO<double>::loadPCD(const string& fileName);

//template
//PointMatcherIO<float>::DataPoints PointMatcherIO<float>::loadPCD(istream& is);
//template
//PointMatcherIO<double>::DataPoints PointMatcherIO<double>::loadPCD(istream& is);

//! @brief Load PCD file
//! @see loadPCD()
template<typename T>
typename PointMatcherIO<T>::DataPoints PointMatcherIO<T>::loadPCD(std::istream& is) {

	//typedef typename DataPoints::Label Label;
	//typedef typename DataPoints::Labels Labels;

	size_t numFields = 0;
	size_t numDataFields = 0; // takes into account the cound of each field for multi row descriptors
	int xFieldCol = -1;
	int yFieldCol = -1;
	int zFieldCol = -1;

	vector<int> descFieldsToKeep;
	map<int,LabelAssociationPair> colToDescPair;
	map<string,int> descLabelToNumRows;
	map<string,int> descLabelToStartingRows;
	vector<int> descDimensions;

	string xFieldType;
	string yFieldType;
	string zFieldType;

	size_t width = 0;
	size_t height = 0;
	size_t numPoints;
	size_t numPointsR = 0; // redundant value specified in POINTS field

	size_t lineNum = 0;

	while (!is.eof())
	{
		string line;
		getline(is, line);

		// get rid of white spaces before/after
		boost::trim (line);

		// ignore comments
		if (line.substr(0,1) == "#")
		{
			lineNum++;
			continue;
		}

		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of("\t\r "), boost::token_compress_on);

		string pcd_version_str;
		if (tokens[0] == "VERSION")
		{
			if (tokens[1] != "0.7" && tokens[1] != ".7")
				throw runtime_error("PCD Parse Error: Only PCD Version 0.7 is supported");
		}

		else if (tokens[0] == "FIELDS")
		{
			numFields = tokens.size() - 1;
			numDataFields = numFields; // in case COUNT is not defined in which case we assume 1 data field per field
			for (size_t i = 1; i < tokens.size(); i++)
			{
				if (tokens[i] == "x")
					xFieldCol = i - 1;
				else if (tokens[i] == "y")
					yFieldCol = i - 1;
				else if (tokens[i] == "z")
					zFieldCol = i - 1;

				else if(descSublabelRegistered(tokens[i]))
				{
					descFieldsToKeep.push_back(i);
					LabelAssociationPair associationPair = getDescAssociationPair(tokens[i]);

					colToDescPair[i] = associationPair;
					descLabelToNumRows[associationPair.second]++;
				}
			}
		}

		else if (tokens[0] == "SIZE")
		{
			if (xFieldCol == -1 || yFieldCol == -1)
				throw runtime_error("PCD Parse Error: x field or y field not defined");
			if (tokens.size()  - 1 !=  numFields)
				throw runtime_error("PCD Parse Error: size not defined for all fields");

//			try {
//				xFieldBytes = boost::lexical_cast<int>(tokens[xFieldCol + 1]);
//				yFieldBytes = boost::lexical_cast<int>(tokens[yFieldCol + 1]);
//				if (zFieldCol > -1)
//					zFieldBytes = boost::lexical_cast<int>(tokens[zFieldCol + 1]);
//			}
//			catch (boost::bad_lexical_cast& e)
//			{
//				throw runtime_error("PCD Parse Error: invalid size field");
//			}

		}

		else if (tokens[0] == "TYPE")
		{
			if (xFieldCol == -1 || yFieldCol == -1)
				throw runtime_error("PCD Parse Error: x field or y field not defined");
			if (tokens.size()  - 1 !=  numFields)
				throw runtime_error("PCD Parse Error: type not defined for all fields");
			xFieldType = tokens[xFieldCol + 1];
			yFieldType = tokens[yFieldCol + 1];

			if (xFieldType != "I" && xFieldType != "U" && xFieldType != "F" &&
					yFieldType != "I" && yFieldType != "U" && yFieldType != "F")
				throw runtime_error("PCD Parse Error: invalid type");

			if (zFieldCol > -1)
			{
				zFieldType = tokens[zFieldCol + 1];
				if (zFieldType != "I" && zFieldType != "U" && zFieldType != "F")
					throw runtime_error("PCD Parse Error: invalid type");
			}
		}

		// overwrite descriptor dimension count with values from header
		else if (tokens[0] == "COUNT")
		{
			if (tokens.size() - 1 != numFields)
				throw runtime_error("PCD Parse Error: COUNT number does not match number of fields");

			// first get total count including fields we aren't using
			numDataFields = 0;

			// we need to overwrite the col to desc pair since there will be more
			// columns now that we have several data counts per field
			map<int, LabelAssociationPair> colToDescPair_ = colToDescPair;
			colToDescPair.clear();


			vector<int>::const_iterator nextFieldToKeepIt = descFieldsToKeep.begin();

			for (size_t i = 1; i < tokens.size(); i++)
			{
				int count = boost::lexical_cast<int>(tokens[i]);
				
				if(descFieldsToKeep.size() != 0)
				{
					if ((int)i == *nextFieldToKeepIt)
					{
						assert(colToDescPair_.find(i) != colToDescPair_.end());

						string descLabel = colToDescPair_[i].second;
						descLabelToNumRows[descLabel] = count;

						for (int p = 0; p < count; p++)
							colToDescPair[numDataFields + p] = LabelAssociationPair(p, descLabel);

						if (nextFieldToKeepIt != descFieldsToKeep.end())
							nextFieldToKeepIt++;
					}
				}

				numDataFields += count;

			}
		}

		else if (tokens[0] == "WIDTH")
		{
			try
			{
				width = boost::lexical_cast<int>(tokens[1]);
			} catch (boost::bad_lexical_cast& e)
			{
				throw runtime_error("PCD Parse Error: invalid width");
			}
		}

		else if (tokens[0] == "HEIGHT")
		{
			try
			{
				height = boost::lexical_cast<int>(tokens[1]);
			} catch (boost::bad_lexical_cast& e)
			{
				throw runtime_error("PCD Parse Error: invalid width");
			}
		}

		// ignore viewpoint for now
		else if (tokens[0] == "VIEWPOINT")
		{
			continue;
		}

		else if (tokens[0] == "POINTS")
		{
			try
			{
				numPointsR = boost::lexical_cast<int>(tokens[1]);
			}
			catch (boost::bad_lexical_cast& e)
			{
				throw runtime_error("PCD Parse Error: invalid number of points");
			}
		}

		else if (tokens[0] == "DATA")
		{
			if (tokens[1] != "ascii")
				throw runtime_error("PCD Parse Error: only ascii data is supported");

			break;
		}

		lineNum++;
	}

	// get number of points
	numPoints = width * height;

	if (numPoints != numPointsR)
		throw runtime_error("PCD Parse Error: POINTS field does not match WIDTH and HEIGHT fields");

	// prepare features matrix
	Matrix features;
	if (zFieldCol > -1)
		features = Matrix(4,numPoints);
	else
		features = Matrix(3,numPoints);

	// Prepare descriptors
	// Do cumulative sum over number of descriptor rows per decriptor to get the starting
	// index row of reach descriptor
	int cumSum = 0;
	for(map<string,int>::const_iterator it = descLabelToNumRows.begin(); it != descLabelToNumRows.end(); it++)
	{
		descLabelToStartingRows[it->first] = cumSum;
		cumSum += it->second;
	}

	// allocate descriptor vectors
	size_t numDescCols = cumSum; // number of descriptor vectors
	Matrix descriptors(numDescCols,numPoints);


	// Now read in the data
	size_t p = 0; // point count
	while (!is.eof())
	{
		string line;
		getline(is, line);

		// get rid of white spaces before/after
		boost::trim (line);

		// ignore comments
		if (line.substr(0,1) == "#")
		{
			lineNum++;
			continue;
		}

		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of("\t\r "), boost::token_compress_on);

		if (tokens.size() != numDataFields)
			throw runtime_error(string("PCD Parse Error: number of data columns does not match number of fields at line: ") + boost::lexical_cast<string>(lineNum));

		features(0,p) = boost::lexical_cast<T>(tokens[xFieldCol]);
		features(1,p) = boost::lexical_cast<T>(tokens[yFieldCol]);

		if (zFieldCol > -1)
		{
			features(2,p) = boost::lexical_cast<float>(tokens[zFieldCol]);
			features(3,p) = 1;
		} else
			features(2,p) = 1;

		for (map<int,LabelAssociationPair>::const_iterator cit = colToDescPair.begin();
				cit != colToDescPair.end(); cit++)
		{
			int startingRow = descLabelToStartingRows[cit->second.second];
			descriptors(startingRow + cit->second.first,p) = boost::lexical_cast<T>(tokens[cit->first]);
		}

		p++;
		lineNum++;

		if (p == numPoints)
			break;

	}

	if (p != numPoints)
	{
		boost::format errorFmt("PCD Parse Error: the number of points in the data %1 is less than the specified number of points %2");
		errorFmt % p % numPoints;
		throw runtime_error(errorFmt.str());
	}

	Labels featureLabels;
	featureLabels.push_back(Label("x"));
	featureLabels.push_back(Label("y"));

	Labels descriptorLabels;
	int n = 0;
	for (map<string,int>::const_iterator it = descLabelToNumRows.begin(); it != descLabelToNumRows.end(); it++)
	{
		descriptorLabels.push_back(Label(it->first,it->second));
		n++;
	}

	if (zFieldCol > -1)
		featureLabels.push_back(Label("z"));

	DataPoints out;

	if (numDescCols > 0)
		out = DataPoints(features, featureLabels, descriptors, descriptorLabels);
	else
		out = DataPoints(features, featureLabels);
	return out;
}

template<typename T>
void PointMatcherIO<T>::savePCD(const DataPoints& data,
		const std::string& fileName) {
	ofstream ofs(fileName.c_str());
	if (!ofs.good())
		throw runtime_error(string("Cannot open file ") + fileName);

	const int pointCount(data.features.cols());
	const int featCount(data.features.rows());
	const int descRows(data.descriptors.rows());
	const int descCount(data.descriptorLabels.size());

	if (pointCount == 0)
	{
		cerr << "Warning, no points, doing nothing" << endl;
		return;
	}

	ofs << "# .PCD v.7 - Point Cloud Data file format\n" <<"VERSION .7\n";
	ofs << "FIELDS";

	for (int f=0; f < (featCount - 1); f++)
	{
		ofs << " " << data.featureLabels[f].text;
	}

	if (descRows == 0)
		ofs << "\n";
	else
	{
		for (int i = 0; i < descCount; i++)
		{
			ofs << " " << data.descriptorLabels[i].text;
		}
		ofs << "\n";
	}

	ofs << "SIZE";
	for (int i =0; i < featCount - 1 + descCount; i++)
	{
		ofs << " 4"; // for float
	}
	ofs << "\n";

	ofs << "TYPE";
	for (int i =0; i < featCount - 1 + descCount; i++)
	{
		ofs << " F"; // for float
	}
	ofs << "\n";

	ofs << "COUNT";
	for (int f = 0; f < featCount - 1 ; f++ )
		ofs << " 1";
	if (descCount == 0)
		ofs << "\n";
	else
	{
		for (int i = 0; i < descCount; i++)
		{
			ofs << " " << data.descriptorLabels[i].span;
		}
		ofs << "\n";
	}

	ofs << "WIDTH " << pointCount << "\n";
	ofs << "HEIGHT 1\n";
	ofs << "POINTS " << pointCount << "\n";
	ofs << "DATA ascii\n";

	// write points
	for (int p = 0; p < pointCount; ++p)
	{
		for (int f = 0; f < featCount - 1; ++f)
		{
			ofs << data.features(f, p);
			if(!(f == featCount-2 && descRows == 0))
				ofs << " ";
		}
		for (int d = 0; d < descRows; ++d)
		{
			ofs << data.descriptors(d, p);
			if(d != descRows-1)
				ofs << " ";
		}
		ofs << "\n";
	}

	ofs.close();
}

template
void PointMatcherIO<float>::savePCD(const DataPoints& data, const std::string& fileName);
template
void PointMatcherIO<double>::savePCD(const DataPoints& data, const std::string& fileName);

template<typename T>
istream & PointMatcherIO<T>::safeGetLine( istream& is, string & t)
{
   t.clear();

       // The characters in the stream are read one-by-one using a std::streambuf.
       // That is faster than reading them one-by-one using the std::istream.
       // Code that uses streambuf this way must be guarded by a sentry object.
       // The sentry object performs various tasks,
       // such as thread synchronization and updating the stream state.

       std::istream::sentry se(is, true);
       std::streambuf* sb = is.rdbuf();

       for(;;) {
           int c = sb->sbumpc();
           switch (c) {
           case '\n':
               return is;
           case '\r':
               if(sb->sgetc() == '\n')
                   sb->sbumpc();
               return is;
           case EOF:
               // Also handle the case when the last line has no line ending
               if(t.empty())
                   is.setstate(std::ios::eofbit);
               return is;
           default:
               t += (char)c;
           }
       }
}


