#pragma once

#ifndef MANET_PROJECT_VISUALIZATION_HEADER_DEF
#define MANET_PROJECT_VISUALIZATION_HEADER_DEF

#include <iostream>
#include <png.h>

#include "ManetDefs.h"

using namespace std;

#define SUCCESS 0

class PngWorker
{
private:
	static png_structp open_png_ptr;
	static png_infop open_info_ptr;
	static FILE* open_img_file;

public:
	static void user_error_fn(png_structp, png_const_charp)
	{
		cout << "[png_create_write_struct]: user_error_fn() called" << endl;
	}

	static void user_warning_fn(png_structp, png_const_charp)
	{
		cout << "[png_create_write_struct]: user_warning_fn() called" << endl;
	}

	static void write_row_callback(png_structp png_ptr, png_uint_32 row, int pass)
	{
		//cout << "write_row_callback() has been called" << endl;
	}

	static int Read(const char * filename, void *** imageData, unsigned int * width, unsigned int * height)
	{
		FILE* imgFile = fopen(filename, "rb");
		if (imgFile == NULL)
		{
			return FOPEN_FAILED;
		}

		png_byte buffer[256];

		fread(buffer, 1, 8, imgFile);
		
		if (png_sig_cmp(buffer, 0, 8))
		{
			return FORMAT_ERROR_NOT_PNG;
		}

		char user_error_ptr[256];

		png_structp png_ptr = png_create_read_struct(
			PNG_LIBPNG_VER_STRING,
			(png_voidp)user_error_ptr,
			user_error_fn,
			user_warning_fn
		);

		if (png_ptr == NULL)
		{
			return PNG_INITIALIZATION_FAILED;
		}

		png_infop info_ptr = png_create_info_struct(png_ptr);

		if (info_ptr == NULL)
		{
			png_destroy_read_struct(&png_ptr, (png_infopp) NULL, (png_infopp) NULL);
			return PNG_INITIALIZATION_FAILED;
		}

		png_infop end_info = png_create_info_struct(png_ptr);

		if (end_info == NULL)
		{
			png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
			return PNG_INITIALIZATION_FAILED;
		}

		if (setjmp(png_jmpbuf(png_ptr)))
		{
			png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
			fclose(imgFile);
			return PNG_READING_ERROR;
		}

		png_init_io(png_ptr, imgFile);
		png_set_sig_bytes(png_ptr, 8);

		png_read_info(png_ptr, info_ptr);

		*width = png_get_image_width(png_ptr, info_ptr);
		*height = png_get_image_height(png_ptr, info_ptr);
		int bit_depth = png_get_bit_depth(png_ptr, info_ptr);
		int color_type = png_get_color_type(png_ptr, info_ptr);
		int filter_method = png_get_filter_type(png_ptr, info_ptr);
		int compression_type = png_get_compression_type(png_ptr, info_ptr);
		int interlace_type = png_get_interlace_type(png_ptr, info_ptr);

		if (!(color_type & PNG_COLOR_MASK_ALPHA))
		{
			cout << "Templates must contain an alpha channel" << endl;
			return FORMAT_ERROR_NO_ALPHA;
		}

		if (bit_depth == 16)
		{
			png_set_strip_16(png_ptr);
		}

		int channels = channels = png_get_channels(png_ptr, info_ptr);
		int rowbytes = png_get_rowbytes(png_ptr, info_ptr);

		*imageData = (void**) new png_bytep[(*height)];
		for (int i = 0; i < (*height); i++)
		{
			(*imageData)[i] = new png_byte[((*width) * 4)]; // ((*width) * 3 == rowbytes)
		}

		png_read_image(png_ptr, (png_bytep*)(*imageData));

		png_read_end(png_ptr, end_info);
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);

		fclose(imgFile);

		return SUCCESS;
	}

	static int WriteEnd()
	{
		png_write_end(open_png_ptr, open_info_ptr);

		png_destroy_write_struct(&open_png_ptr, &open_info_ptr);

		fclose(open_img_file);

		return SUCCESS;
	}

	static int WriteRow(unsigned char * data)
	{
		png_write_row(open_png_ptr, (png_bytep) data);

		return SUCCESS;
	}

	static int WriteOpen(const char * fileName, int width, int height)
	{
		open_img_file = fopen(fileName, "wb");

		if (open_img_file == NULL)
		{
			cout << "[PngWriter.Test()]: fopen() failed";
			return 1;
		}

		char user_error_ptr[280];

		open_png_ptr = png_create_write_struct(
			PNG_LIBPNG_VER_STRING,
			(png_voidp)user_error_ptr,
			user_error_fn,
			user_warning_fn
		);

		if (!open_png_ptr)
		{
			return PNG_INITIALIZATION_FAILED;
		}

		open_info_ptr = png_create_info_struct(open_png_ptr);

		if (!open_info_ptr)
		{
			png_destroy_write_struct(&open_png_ptr,
				(png_infopp)NULL);
			return PNG_INITIALIZATION_FAILED;
		}

		if (setjmp(png_jmpbuf(open_png_ptr)))
		{
			png_destroy_write_struct(&open_png_ptr, &open_info_ptr);
			fclose(open_img_file);
			return PNG_WRITING_ERROR;
		}

		png_init_io(open_png_ptr, open_img_file);

		png_set_write_status_fn(open_png_ptr, write_row_callback);

		png_set_IHDR(
			open_png_ptr,
			open_info_ptr,
			(png_uint_32) width, // width
			(png_uint_32) height, // height
			8, // bit depth
			PNG_COLOR_TYPE_RGB_ALPHA,
			PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_DEFAULT,
			PNG_FILTER_TYPE_DEFAULT
		);

		png_write_info(open_png_ptr, open_info_ptr);
		// png_set_filler(png_ptr, 0, PNG_FILLER_BEFORE);
		// png_set_shift(png_ptr, &sig_bit);
		// if(bit_depth > 8) png_set_swap(png_ptr);
		// if(bit_depth < 8) png_set_packswap(png_ptr);
		//png_set_bgr(png_ptr);
		//png_set_flush(png_ptr, 10);

		return 0;
	}
};

png_structp PngWorker::open_png_ptr;
png_infop PngWorker::open_info_ptr;
FILE* PngWorker::open_img_file;

class GraphVisualizer
{
public:
	/*
	 * A struct to be passed into DrawGraphRectangularGrid
	 * 
	 * const char ** digits - holds an array of 10 const char *
	 *     each of which stores a filename of the PNG image
	 *     of the respective (with respect to an array index) digit
	 * 
	 */
	struct RectangularGridTemplates
	{
		vector<string> digits;
		string gridCross;
		string gridCrossX0;
		string gridCrossY0;
		string gridCrossZero;
		string gridArrowXAxis;
		string gridArrowYAxis;
		string gridEmpty;

		int playersCount;

		vector<string> staticAgents;
		vector<string> drones;

		RectangularGridTemplates(int playersCount) : playersCount(playersCount)
		{
			digits.resize(10);
			staticAgents.resize(playersCount);
			drones.resize(playersCount);
		}

		RectangularGridTemplates() : playersCount(0)
		{
			digits.resize(10);
		}

		~RectangularGridTemplates() {}
	};

	struct ImageData
	{
		void ** data;
		uint32 width;
		uint32 height;

		ImageData() {}
	};

	class PosComparator
	{
	public:
		bool operator()(pair<int, int> lhs, pair<int, int> rhs)
		{
			return (lhs.second < rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first));
		}
	};

	static bool FindPos(vector<pair<int,int>> &vec, pair<int, int> toFind)
	{
		PosComparator pc;
		return binary_search(vec.begin(), vec.end(), toFind, pc);
	}

	/*
	* this function assumes that dst alpha equals 1
	*/
	static void AddDataWithAlpha(byte * dst, byte * src, uint32 lengthPixels)
	{
		for (uint32 i = 0; i < lengthPixels; i++)
		{
			uint16 alphaSrc = (uint16)src[i * 4 + 3];
			uint16 rSrc = (uint16)src[i * 4] * alphaSrc;
			uint16 gSrc = (uint16)src[i * 4 + 1] * alphaSrc;
			uint16 bSrc = (uint16)src[i * 4 + 2] * alphaSrc;

			uint16 alphaDst = (uint16)(0xFF - src[i * 4 + 3]);
			uint16 rDst = (uint16)dst[i * 4] * alphaDst;
			uint16 gDst = (uint16)dst[i * 4 + 1] * alphaDst;
			uint16 bDst = (uint16)dst[i * 4 + 2] * alphaDst;

			dst[i * 4] = (byte)((rSrc + rDst) / (uint16)0xFF);
			dst[i * 4 + 1] = (byte)((gSrc + gDst) / (uint16)0xFF);
			dst[i * 4 + 2] = (byte)((bSrc + bDst) / (uint16)0xFF);
			dst[i * 4 + 3] = 0xFF;
		}
	}

	static void DrawGraphRectangularGrid(
		const char* outFileName,
		RectangularGridTemplates& templates,
		int playersCount,
		vector<vector<pair<int, int>>>& staticAgentsPositions,
		vector<vector<pair<int, int>>>& dronesPositions
	)
	{
		ImageData* digits = new ImageData[10];
		ImageData gridCross;
		ImageData gridCrossX0;
		ImageData gridCrossY0;
		ImageData gridCrossZero;
		ImageData gridArrowXAxis;
		ImageData gridArrowYAxis;
		ImageData gridEmpty;

		ImageData* staticAgents = new ImageData[playersCount];
		ImageData* drones = new ImageData[playersCount];

		for (int i = 0; i < 10; i++)
		{
			PngWorker::Read(templates.digits[0].c_str(), &digits[0].data, &digits[0].width, &digits[0].height);
		}

		PngWorker::Read(templates.gridArrowYAxis.c_str(), &gridArrowYAxis.data, &gridArrowYAxis.width, &gridArrowYAxis.height);
		PngWorker::Read(templates.gridCross.c_str(), &gridCross.data, &gridCross.width, &gridCross.height);
		PngWorker::Read(templates.gridCrossX0.c_str(), &gridCrossX0.data, &gridCrossX0.width, &gridCrossX0.height);
		PngWorker::Read(templates.gridCrossY0.c_str(), &gridCrossY0.data, &gridCrossY0.width, &gridCrossY0.height);
		PngWorker::Read(templates.gridCrossZero.c_str(), &gridCrossZero.data, &gridCrossZero.width, &gridCrossZero.height);
		PngWorker::Read(templates.gridArrowXAxis.c_str(), &gridArrowXAxis.data, &gridArrowXAxis.width, &gridArrowXAxis.height);
		PngWorker::Read(templates.gridEmpty.c_str(), &gridEmpty.data, &gridEmpty.width, &gridEmpty.height);

		for (int i = 0; i < playersCount; i++)
		{
			PngWorker::Read(templates.staticAgents[i].c_str(), &staticAgents[i].data, &staticAgents[i].width, &staticAgents[i].height);
			PngWorker::Read(templates.drones[i].c_str(), &drones[i].data, &drones[i].width, &drones[i].height);
		}

		PosComparator pc;

		for (int i = 0; i < playersCount; i++)
		{
			sort(staticAgentsPositions[i].begin(), staticAgentsPositions[i].end(), pc);
			sort(dronesPositions[i].begin(), dronesPositions[i].end(), pc);
		}

		int top = INT_MIN;
		int left = INT_MAX;
		int bottom = INT_MAX;
		int right = INT_MIN;

		for (int i = 0; i < playersCount; i++)
		{
			for (int j = 0; j < staticAgentsPositions[i].size(); j++)
			{
				top = max(top, staticAgentsPositions[i][j].second);
				bottom = min(bottom, staticAgentsPositions[i][j].second);
				right = max(right, staticAgentsPositions[i][j].first);
				left = min(left, staticAgentsPositions[i][j].first);
			}

			for (int j = 0; j < dronesPositions[i].size(); j++)
			{
				top = max(top, dronesPositions[i][j].second);
				bottom = min(bottom, dronesPositions[i][j].second);
				right = max(right, dronesPositions[i][j].first);
				left = min(left, dronesPositions[i][j].first);
			}
		}

		/*
		values below show how many empty rows and columns
		are added between axes and bottom left corner of
		the "graph rectangle" (rectangle of top left - bottom right)
		vertices
		*/
		//unsigned int bottom_auto_spacing = 

		unsigned int width = (right + 2) * gridCross.width;
		unsigned int height = (top + 2) * gridCross.height;

		PngWorker::WriteOpen(outFileName, width, height);

		unsigned char* row = new unsigned char[width*4];

		for (int row_number = 0; row_number < height; row_number++)
		{
			for (int x = 0; x < (right + 2); x++)
			{
				int32 y = (height - row_number - 1) / gridCross.height;
				int pattern_row_number = row_number % gridCross.height;

				byte ** patternData = (byte**) gridCross.data;

				uint32 cellWidthPixels = gridCross.width;
				uint32 cellWidthBytes = cellWidthPixels * 4;

				if (x == 0 && y == 0)
				{
					patternData = (byte**) gridCrossZero.data;
				}
				else
				{
					if (x == 0)
					{
						if (y == (top + 1)) { patternData = (byte**) gridArrowYAxis.data; }
						else { patternData = (byte**) gridCrossX0.data; }
					}
					else
					{
						if (y == (top + 1)) { patternData = (byte**) gridEmpty.data; }
					}

					if (y == 0)
					{
						if (x == right + 1) { patternData = (byte**) gridArrowXAxis.data; }
						else { patternData = (byte**) gridCrossY0.data; }
					}
					else
					{
						if (x == right + 1) { patternData = (byte**) gridEmpty.data; }
					}
				}

				memcpy(row + (x * cellWidthBytes), patternData[pattern_row_number], cellWidthBytes);

				for (int playerNumber = 0; playerNumber < playersCount; playerNumber++)
				{
					if (FindPos(staticAgentsPositions[playerNumber], make_pair(x, y)))
					{
						patternData = (byte**) staticAgents[playerNumber].data;
						AddDataWithAlpha(row + (x * cellWidthBytes), patternData[pattern_row_number], cellWidthPixels);
					}
				}

				for (int playerNumber = 0; playerNumber < playersCount; playerNumber++)
				{
					if (FindPos(dronesPositions[playerNumber], make_pair(x, y)))
					{
						patternData = (byte**) drones[playerNumber].data;
						AddDataWithAlpha(row + (x * cellWidthBytes), patternData[pattern_row_number], cellWidthPixels);
					}
				}
			}

			PngWorker::WriteRow(row);
		}

		PngWorker::WriteEnd();
	}
};

#endif /* MANET_PROJECT_VISUALIZATION_HEADER_DEF */
