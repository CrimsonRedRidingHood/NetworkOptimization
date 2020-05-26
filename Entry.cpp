#include "Network.hpp"

#include <map>
#include <set>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

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
	bool operator()(const Pos& lhs, const Pos& rhs)
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
	bool operator()(const Agent& lhs, const Agent& rhs)
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
	int operator()(const Agent& lhs, const Agent& rhs)
	{
		int idist = abs(lhs.position.i - rhs.position.i);
		int jdist = abs(lhs.position.j - rhs.position.j);
		return (idist + jdist);
	}
};

int main()
{
	Network<Pos, Agent, DoNeighbor, int, CalculateDistance, PosComp> network;
	
	const int player1AgentsNumber = 8;
	const int player2AgentsNumber = 8;

	int player1Pos[player1AgentsNumber][2] = {
		/*
		{1,1},{2,1},{2,2},
		{2,3},{2,4},{3,4},
		{4,4},{4,3},{4,2},{5,2}
		*/
		{2,1},{2,2},{3,2},
		{3,3},{4,2},{4,3},
		{5,3},{6,3}
	};

	int player2Pos[player2AgentsNumber][2] = {
		/*
		{3,5},{4,5},{5,5},
		{6,5},{4,4},{4,3},
		{4,2},{5,2},{6,2},
		{7,2},{8,2}
		*/
		{1,4},{1,5},{2,5},
		{3,5},{3,6},{4,5},
		{4,4},{4,3}
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

	cout << "(" << solution.first.i << ", " << solution.first.j << ") (" <<
		solution.second.i << ", " << solution.second.j << ")" << endl;

	cout << "Press Enter to quit...";

	fgetc( stdin );
	
	return 0;
}