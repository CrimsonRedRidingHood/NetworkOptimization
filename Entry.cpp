#include "Network.hpp"
#include "Visualization.hpp"

#include <map>
#include <set>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "ManetDefs.h"

using namespace std;

class Pos
{
public:
	int i;
	int j;

	Pos() : i(0), j(0) {}

	Pos(int iInitial, int jInitial) : i(iInitial), j(jInitial) { }

	Pos& operator=(Pos rhs)
	{
		i = rhs.i;
		j = rhs.j;
		return *this;
	}

	Pos(const Pos& another) : i(another.i), j(another.j) {}

	Pos Left() { return Pos(this->i - 1, this->j); }
	Pos Right() { return Pos(this->i + 1, this->j); }
	Pos Up() { return Pos(this->i, this->j - 1); }
	Pos Down() { return Pos(this->i, this->j + 1); }

	bool operator<(const Pos another) const
	{
		return(this->i < another.i) || ((this->i == another.i) && (this->j < another.j));
	}
};

class PosComp
{
public:
	bool operator()(const Pos& lhs, const Pos& rhs) const
	{
		return (lhs.i < rhs.i) || ((lhs.i == rhs.i) && (lhs.j < rhs.j));
	}
};

class Agent
{
public:
	Pos position;

	int playerID;

	~Agent() {}

	Agent() : position(), playerID(-1) {}

	Agent(int iInitial, int jInitial, int playerID) : position(iInitial, jInitial), playerID(playerID) {}

	Agent(Pos positionInitial, int playerID) : position(positionInitial), playerID(playerID) {}

	Agent(int playerID) : playerID(playerID) {}

	Agent(Agent& another) : position(another.position), playerID(another.playerID) {}
};

class DoNeighbor
{
public:
	bool operator()(const Agent& lhs, const Agent& rhs) const
	{
		int idist = abs(lhs.position.i - rhs.position.i);
		int jdist = abs(lhs.position.j - rhs.position.j);

		bool distanceIsRight = ((idist + jdist) <= 1) && ((idist + jdist) > 0);
		bool eitherIsADrone = ((lhs.playerID | rhs.playerID) & 0x10) != 0;
		bool samePlayer = (lhs.playerID == rhs.playerID);

		//return (distanceIsRight && (eitherIsADrone || samePlayer));
		return distanceIsRight;
	}
};

class CalculateDistance
{
public:
	int operator()(const Agent& lhs, const Agent& rhs) const
	{
		int idist = abs(lhs.position.i - rhs.position.i);
		int jdist = abs(lhs.position.j - rhs.position.j);
		return (idist + jdist);
	}
};

int TestRead(const char* filename, unsigned char*** imageData, unsigned int* width, unsigned int* height)
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
		NULL,
		NULL
	);

	if (png_ptr == NULL)
	{
		return PNG_INITIALIZATION_FAILED;
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);

	if (info_ptr == NULL)
	{
		png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
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

	int channels = channels = png_get_channels(png_ptr, info_ptr);
	int rowbytes = png_get_rowbytes(png_ptr, info_ptr);

	if (color_type & PNG_COLOR_MASK_ALPHA)
	{
		png_set_strip_alpha(png_ptr);
	}

	*imageData = new png_bytep[(*height)];
	for (int i = 0; i < (*height); i++)
	{
		(*imageData)[i] = new png_byte[((*width) * 3)]; // ((*width) * 3 == rowbytes)
	}

	png_read_rows(png_ptr, (png_bytep*)(*imageData), NULL, (*height));

	png_read_end(png_ptr, end_info);
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);

	fclose(imgFile);

	return SUCCESS;
}

int main()
{
	/*
	void ** img1 = NULL;
	unsigned int width1;
	unsigned int height1;

	void ** img2 = NULL;
	unsigned int width2;
	unsigned int height2;

	PngWorker::Read("Visualization\\RectangularGrid\\player_1_static_agent.png", &img1, &width1, &height1);
	PngWorker::Read("Visualization\\RectangularGrid\\player_2_drone.png", &img2, &width2, &height2);

	PngWorker::WriteOpen("test.png", width1 + width2, height1);

	unsigned char* row = new unsigned char[(width1 + width1) * 3];

	for (int i = 0; i < height1; i++)
	{
		memcpy(row, img1[i], width1 * 3);
		memcpy(row + (width1 * 3), img2[i], width2 * 3);
		PngWorker::WriteRow(row);
	}

	PngWorker::WriteEnd();

	for (int i = 0; i < height1; i++)
	{
		delete[] img1[i];
	}

	delete[] img1;

	for (int i = 0; i < height2; i++)
	{
		delete[] img2[i];
	}

	delete[] img2;

	return 0;
	*/

	Network<Pos, Agent, DoNeighbor, int, CalculateDistance, PosComp> network;
	
	const int player1AgentsNumber = 13;
	const int player2AgentsNumber = 17;

	int player1Pos[player1AgentsNumber][2] = {
		/*
		{1,1},{2,1},{2,2},
		{2,3},{2,4},{3,4},
		{4,4},{4,3},{4,2},{5,2}
		*/
		{2,3},{2,4},{2,5},
		{2,6},{2,7},{2,8},
		{2,9},{2,10},{3,10},
		{4,10},{4,9},{4,8},
		{4,7}
	};

	int player2Pos[player2AgentsNumber][2] = {
		/*
		{3,5},{4,5},{5,5},
		{6,5},{4,4},{4,3},
		{4,2},{5,2},{6,2},
		{7,2},{8,2}
		*/
		{5,3},{6,3},{7,3},
		{8,3},{9,3},{10,3},
		{11,3},{11,4},{11,5},
		{10,5},{9,5},{8,5},
		{7,5},{6,5},{5,5},
		{4,5},{3,5}
	};

	Agent ** player1Agents = new Agent*[player1AgentsNumber];
	Agent ** player2Agents = new Agent*[player2AgentsNumber];

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		player1Agents[i] = new Agent(player1Pos[i][0], player1Pos[i][1], 0x01);
		network.Add(player1Agents[i]->position, player1Agents[i]);
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		player2Agents[i] = new Agent(player2Pos[i][0], player2Pos[i][1], 0x02);
		network.Add(player2Agents[i]->position, player2Agents[i]);
	}

	int player1Diameter = network.GetDiameter(
		[](Agent x) {return (x.playerID == 0x01);},
		[](Agent x) {return (x.playerID == 0x01);},
		LONG_MAX / 2
	);

	int player2Diameter = network.GetDiameter(
		[](Agent x) {return (x.playerID == 0x02);},
		[](Agent x) {return (x.playerID == 0x02);},
		LONG_MAX / 2
	);

	set<Pos, PosComp> existingPositions;
	set<Pos, PosComp> candidatesUniquifier;

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		candidatesUniquifier.insert(Pos(player1Pos[i][0] - 1, player1Pos[i][1]));
		candidatesUniquifier.insert(Pos(player1Pos[i][0] + 1, player1Pos[i][1]));
		candidatesUniquifier.insert(Pos(player1Pos[i][0], player1Pos[i][1] - 1));
		candidatesUniquifier.insert(Pos(player1Pos[i][0], player1Pos[i][1] + 1));

		existingPositions.insert(Pos(player1Pos[i][0], player1Pos[i][1]));
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		candidatesUniquifier.insert(Pos(player2Pos[i][0] - 1, player2Pos[i][1]));
		candidatesUniquifier.insert(Pos(player2Pos[i][0] + 1, player2Pos[i][1]));
		candidatesUniquifier.insert(Pos(player2Pos[i][0], player2Pos[i][1] - 1));
		candidatesUniquifier.insert(Pos(player2Pos[i][0], player2Pos[i][1] + 1));

		existingPositions.insert(Pos(player2Pos[i][0], player2Pos[i][1]));
	}

	vector<Pos> candidates;

	for (set<Pos, PosComp>::iterator iter = candidatesUniquifier.begin(); iter != candidatesUniquifier.end(); iter++)
	{
		candidates.push_back(*iter);
	}

	CalculateDistance distanceCalculator;

	Agent * player1Drone = new Agent(0x11);
	Agent * player2Drone = new Agent(0x12);

	map<int, set<pair<Pos, Pos>>> paretoOptimal;

	for (int i = 0; i < candidates.size(); i++)
	{
		for (int j = i; j < candidates.size(); j++)
		{
			/*
				checking if the pair of positions
				can possibly provide some improvement
			*/

			int idist = abs(candidates[i].i - candidates[j].i);
			int jdist = abs(candidates[i].j - candidates[j].j);

			int player1NeighborsNumber =
				(existingPositions.find(candidates[i].Left()) != existingPositions.end()) +
				(existingPositions.find(candidates[i].Right()) != existingPositions.end()) +
				(existingPositions.find(candidates[i].Up()) != existingPositions.end()) + 
				(existingPositions.find(candidates[i].Down()) != existingPositions.end());
			
			int player2NeighborsNumber =
				(existingPositions.find(candidates[j].Left()) != existingPositions.end()) +
				(existingPositions.find(candidates[j].Right()) != existingPositions.end()) +
				(existingPositions.find(candidates[j].Up()) != existingPositions.end()) +
				(existingPositions.find(candidates[j].Down()) != existingPositions.end());

			bool connected = (idist + jdist) == 1;
			bool haveNeighbors = ((player1NeighborsNumber > 0) && (player2NeighborsNumber > 0));
			bool haveMoreThanOneNeighbor = ((player1NeighborsNumber > 1) && (player2NeighborsNumber > 1));

			if ((haveNeighbors && connected) || haveMoreThanOneNeighbor)
			{
				player1Drone->position = candidates[i];
				player2Drone->position = candidates[j];
				network.Add(candidates[i], player1Drone);
				network.Add(candidates[j], player2Drone);

				int player1Gain =
					player1Diameter - network.GetDiameter(
						[](Agent x) {return ((x.playerID == 0x01) || ((x.playerID & 0x10) != 0));},
						[](Agent x) {return (x.playerID == 0x01);},
						LONG_MAX / 2
					);

				int player2Gain =
					player2Diameter - network.GetDiameter(
						[](Agent x) {return ((x.playerID == 0x02) || ((x.playerID & 0x10) != 0));},
						[](Agent x) {return (x.playerID == 0x02);},
						LONG_MAX / 2
					);

				/*
					At this point we have the drones set
					and diameter calculated so we can perform
					any computations needed to figure out if
					this set of positions is valuable or not
				*/

				// finding Pareto optimal solutions
				int newResult = player1Gain + player2Gain;

				if (paretoOptimal.size() != 0)
				{
					if (paretoOptimal.begin()->first <= newResult)
					{
						if (paretoOptimal.begin()->first < newResult)
						{
							paretoOptimal.clear();
						}

						paretoOptimal[newResult].insert(make_pair(player1Drone->position, player2Drone->position));
					}
				}
				else
				{
					paretoOptimal[newResult].insert(make_pair(player1Drone->position, player2Drone->position));
				}

				network.Remove(candidates[i], player1Drone);
				network.Remove(candidates[j], player2Drone);
			}
		}
	}

	double bestPageRank = 0.0;
	pair<Pos, Pos> solution;

	for (map<int, set<pair<Pos, Pos>>>::iterator iter = paretoOptimal.begin(); iter != paretoOptimal.end(); iter++)
	{
		for (set<pair<Pos, Pos>>::iterator setiter = iter->second.begin(); setiter != iter->second.end(); setiter++)
		{
			Pos p1pos = setiter->first;
			Pos p2pos = setiter->second;

			player1Drone->position = p1pos;
			player2Drone->position = p2pos;
			network.Add(p1pos, player1Drone);
			network.Add(p2pos, player2Drone);

			vector<pair<Agent*, double>> pageRank = network.GetPageRank(
				[](Agent) {return true;},
				[](Agent) {return true;}
			);

			double currentPageRank = 0.0;

			for (int i = 0; i < pageRank.size(); i++)
			{
				if (pageRank[i].first == player1Drone || pageRank[i].first == player2Drone)
				{
					currentPageRank += pageRank[i].second;
				}
			}

			if (currentPageRank > bestPageRank)
			{
				bestPageRank = currentPageRank;
				solution = make_pair(p1pos, p2pos);
			}

			network.Remove(p1pos, player1Drone);
			network.Remove(p2pos, player2Drone);
		}
	}

	player1Drone->position = solution.first;
	player2Drone->position = solution.second;

	network.Add(solution.first, player1Drone);
	network.Add(solution.second, player2Drone);

	ofstream outfile("output.txt");

	vector<vector<double>> adjMatrix;
	vector<Agent*> indexedAgents;

	const int numberOfPlayers = 2;

	outfile << numberOfPlayers << endl;

	for (int playerNumber = 1; playerNumber <= numberOfPlayers; playerNumber++)
	{
		network.GenerateAdjacencyMatrix(
			[playerNumber](Agent x) {return ((x.playerID == playerNumber) || ((x.playerID & 0x10) != 0));},
			adjMatrix,
			indexedAgents
		);

		outfile << indexedAgents.size() << endl;

		for (int i = 0; i < indexedAgents.size(); i++)
		{
			outfile << indexedAgents[i]->position.i << " " << indexedAgents[i]->position.j << " " << indexedAgents[i]->playerID << endl;
		}

		for (int i = 0; i < adjMatrix.size(); i++)
		{
			for (int j = 0; j < adjMatrix.size(); j++)
			{
				outfile << adjMatrix[i][j] << " ";
			}

			outfile << endl;
		}

		indexedAgents.clear();

		for (int i = 0; i < adjMatrix.size(); i++)
		{
			adjMatrix[i].clear();
		}
		adjMatrix.clear();
	}

	outfile.close();

	vector<vector<pair<int, int>>> drones;
	drones.resize(2);
	drones[0].push_back(make_pair(player1Drone->position.i, player1Drone->position.j));
	drones[1].push_back(make_pair(player2Drone->position.i, player2Drone->position.j));

	delete player1Drone;
	delete player2Drone;

	vector<vector<pair<int, int>>> staticAgents;
	staticAgents.resize(2);
	for (int i = 0; i < player1AgentsNumber; i++)
	{
		staticAgents[0].push_back(make_pair(player1Agents[i]->position.i, player1Agents[i]->position.j));
	}
	for (int i = 0; i < player2AgentsNumber; i++)
	{
		staticAgents[1].push_back(make_pair(player2Agents[i]->position.i, player2Agents[i]->position.j));
	}

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		network.Remove(player1Agents[i]->position, player1Agents[i]);
		delete player1Agents[i];
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		network.Remove(player2Agents[i]->position, player2Agents[i]);
		delete player2Agents[i];
	}

	delete[] player1Agents;
	delete[] player2Agents;

	GraphVisualizer::RectangularGridTemplates templates(2);

	string imgsFolder = "Visualization\\RectangularGrid\\";

	templates.gridCross = imgsFolder + "grid_cross.png";
	templates.gridCrossX0 = imgsFolder + "grid_cross_x0.png";
	templates.gridCrossY0 = imgsFolder + "grid_cross_y0.png";
	templates.gridCrossZero = imgsFolder + "grid_cross_zero.png";
	templates.gridArrowXAxis = imgsFolder + "grid_arrow_x_axis.png";
	templates.gridArrowYAxis = imgsFolder + "grid_arrow_y_axis.png";
	templates.gridEmpty = imgsFolder + "grid_empty.png";

	templates.playersCount = 2;

	for (int i = 0; i < templates.playersCount; i++)
	{
		templates.staticAgents[i] = imgsFolder + "player_" + ((char)(i + 1 + '0')) + "_static_agent.png";
		templates.drones[i] = imgsFolder + "player_" + ((char)(i + 1 + '0')) + "_drone.png";
	}

	for (int i = 0; i < 10; i++)
	{
		templates.digits[i] = imgsFolder + "digit_" + ((char)(i + '0')) + ".png";
	}

	GraphVisualizer::DrawGraphRectangularGrid(
		"result.png",
		templates,
		2,
		staticAgents,
		drones
	);
	
	return 0;
}