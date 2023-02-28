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

class Pos_RectangularGrid
{
public:
	int i;
	int j;

	Pos_RectangularGrid() : i(0), j(0) {}

	Pos_RectangularGrid(int iInitial, int jInitial) : i(iInitial), j(jInitial) { }

	Pos_RectangularGrid& operator=(Pos_RectangularGrid rhs)
	{
		i = rhs.i;
		j = rhs.j;
		return *this;
	}

	Pos_RectangularGrid(const Pos_RectangularGrid& another) : i(another.i), j(another.j) {}

	Pos_RectangularGrid Left() { return Pos_RectangularGrid(this->i - 1, this->j); }
	Pos_RectangularGrid Right() { return Pos_RectangularGrid(this->i + 1, this->j); }
	Pos_RectangularGrid Up() { return Pos_RectangularGrid(this->i, this->j - 1); }
	Pos_RectangularGrid Down() { return Pos_RectangularGrid(this->i, this->j + 1); }

	bool operator<(const Pos_RectangularGrid another) const
	{
		return(this->i < another.i) || ((this->i == another.i) && (this->j < another.j));
	}
};

class PosComp
{
public:
	bool operator()(const Pos_RectangularGrid& lhs, const Pos_RectangularGrid& rhs) const
	{
		return (lhs.i < rhs.i) || ((lhs.i == rhs.i) && (lhs.j < rhs.j));
	}
};

class Agent
{
public:
	Pos_RectangularGrid position;

	int playerID;

	~Agent() {}

	Agent() : position(), playerID(-1) {}

	Agent(int iInitial, int jInitial, int playerID) : position(iInitial, jInitial), playerID(playerID) {}

	Agent(Pos_RectangularGrid positionInitial, int playerID) : position(positionInitial), playerID(playerID) {}

	Agent(int playerID) : playerID(playerID) {}

	Agent(Agent& another) : position(another.position), playerID(another.playerID) {}
};

class DoNeighbor_RectangularGrid
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

class CalculateDistance_Manhattan
{
public:
	int operator()(const Agent& lhs, const Agent& rhs) const
	{
		int idist = abs(lhs.position.i - rhs.position.i);
		int jdist = abs(lhs.position.j - rhs.position.j);
		return (idist + jdist);
	}
};

void buildCandidatesAndExistingSets_RectangularGrid(
	vector<vector<pair<int,int>>> players,
	vector<Pos_RectangularGrid>& candidates,
	set<Pos_RectangularGrid, PosComp>& existingPositions
)
{
	set<Pos_RectangularGrid, PosComp> candidatesUniquifier;

	candidates.clear();
	existingPositions.clear();

	for (int i = 0; i < players.size(); i++)
	{
		for (int j = 0; j < players[i].size(); j++)
		{
			Pos_RectangularGrid current = Pos_RectangularGrid(players[i][j].first, players[i][j].second);

			candidatesUniquifier.insert(current.Left());
			candidatesUniquifier.insert(current.Right());
			candidatesUniquifier.insert(current.Up());
			candidatesUniquifier.insert(current.Down());

			existingPositions.insert(current);
		}
	}

	for (
		set<Pos_RectangularGrid, PosComp>::iterator iter = candidatesUniquifier.begin();
		iter != candidatesUniquifier.end();
		iter++
		)
	{
		candidates.push_back(*iter);
	}
}

template<typename... NetworkTypes>
void writeResultsToFile(string filename, Network<NetworkTypes...>& network)
{
	ofstream outfile("output.txt");

	vector<vector<double>> adjMatrix;
	vector<Agent*> indexedAgents;

	const int numberOfPlayers = 2;

	outfile << numberOfPlayers << endl;

	for (int playerNumber = 1; playerNumber <= numberOfPlayers; playerNumber++)
	{
		network.GenerateAdjacencyMatrix(
			[playerNumber](Agent x) {return ((x.playerID == playerNumber) || ((x.playerID & 0x10) != 0)); },
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
}

void visualizeResults_RectangularGrid(
	const char* outFilename,
	int playersCount,
	vector<vector<pair<int, int>>>& staticAgents,
	vector<vector<pair<int, int>>>& drones
)
{
	GraphVisualizer::RectangularGridTemplates templates(playersCount);

	string imgsFolder = "Visualization\\RectangularGrid\\";

	templates.gridCross = imgsFolder + "grid_cross.png";
	templates.gridCrossX0 = imgsFolder + "grid_cross_x0.png";
	templates.gridCrossY0 = imgsFolder + "grid_cross_y0.png";
	templates.gridCrossZero = imgsFolder + "grid_cross_zero.png";
	templates.gridArrowXAxis = imgsFolder + "grid_arrow_x_axis.png";
	templates.gridArrowYAxis = imgsFolder + "grid_arrow_y_axis.png";
	templates.gridEmpty = imgsFolder + "grid_empty.png";

	templates.playersCount = playersCount;

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
		outFilename,
		templates,
		playersCount,
		staticAgents,
		drones
	);
}

int main()
{
	Network<Pos_RectangularGrid, Agent, DoNeighbor_RectangularGrid, int, CalculateDistance_Manhattan, PosComp> network;
	
	int player1AgentsNumber = 13;
	int player2AgentsNumber = 17;

	vector<pair<int, int>> player1Pos;
	vector<pair<int, int>> player2Pos;
	
	ifstream inpfile;
	inpfile.open("input.txt");

	/*
	int player1Pos[player1AgentsNumber][2] = {
		/old/
		{1,1},{2,1},{2,2},
		{2,3},{2,4},{3,4},
		{4,4},{4,3},{4,2},{5,2}
		/old-end/
		{2,3},{2,4},{2,5},
		{2,6},{2,7},{2,8},
		{2,9},{2,10},{3,10},
		{4,10},{4,9},{4,8},
		{4,7}
	};

	int player2Pos[player2AgentsNumber][2] = {
		/old/
		{3,5},{4,5},{5,5},
		{6,5},{4,4},{4,3},
		{4,2},{5,2},{6,2},
		{7,2},{8,2}
		/old-end/
		{5,3},{6,3},{7,3},
		{8,3},{9,3},{10,3},
		{11,3},{11,4},{11,5},
		{10,5},{9,5},{8,5},
		{7,5},{6,5},{5,5},
		{4,5},{3,5}
	};
	*/

	inpfile >> player1AgentsNumber;

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		int x, y;
		inpfile >> x >> y;
		player1Pos.push_back(make_pair(x,y));
	}

	inpfile >> player2AgentsNumber;

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		int x, y;
		inpfile >> x >> y;
		player2Pos.push_back(make_pair(x, y));
	}

	inpfile.close();

	Agent ** player1Agents = new Agent*[player1AgentsNumber];
	Agent ** player2Agents = new Agent*[player2AgentsNumber];

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		player1Agents[i] = new Agent(player1Pos[i].first, player1Pos[i].second, 0x01);
		network.Add(player1Agents[i]->position, player1Agents[i]);
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		player2Agents[i] = new Agent(player2Pos[i].first, player2Pos[i].second, 0x02);
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

	set<Pos_RectangularGrid, PosComp> existingPositions;
	vector<vector<pair<int, int>>> staticAgents;
	vector<Pos_RectangularGrid> candidates;

	// add players' positions to "players" vector in order to feed them to the function below
	staticAgents.push_back(player1Pos);
	staticAgents.push_back(player2Pos);

	buildCandidatesAndExistingSets_RectangularGrid(staticAgents, candidates, existingPositions);

	Agent * player1Drone = new Agent(0x11);
	Agent * player2Drone = new Agent(0x12);

	/*
	* the map below holds the pair of values:
	* (best value (same for all elements), solution that yields the value)
	*/
	map<int, set<pair<Pos_RectangularGrid, Pos_RectangularGrid>>> paretoOptimal;

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

	/*
	* Here, we seek for the best solution
	* with respect to the PageRank measure
	*/
	double bestPageRank = 0.0;
	pair<Pos_RectangularGrid, Pos_RectangularGrid> solution;

	for (auto iter = paretoOptimal.begin(); iter != paretoOptimal.end(); iter++)
	{
		for (auto setiter = iter->second.begin(); setiter != iter->second.end(); setiter++)
		{
			Pos_RectangularGrid p1pos = setiter->first;
			Pos_RectangularGrid p2pos = setiter->second;

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

	/*
	* Writing the solution in the corresponding file
	*/
	writeResultsToFile("output.txt", network);

	vector<vector<pair<int, int>>> drones;
	drones.resize(2);
	drones[0].push_back(make_pair(player1Drone->position.i, player1Drone->position.j));
	drones[1].push_back(make_pair(player2Drone->position.i, player2Drone->position.j));

	// cleanup phase
	network.Remove(solution.first, player1Drone);
	network.Remove(solution.second, player2Drone);

	delete player1Drone;
	delete player2Drone;

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

	// visualization
	visualizeResults_RectangularGrid("result.png", 2, staticAgents, drones);
	
	return 0;
}