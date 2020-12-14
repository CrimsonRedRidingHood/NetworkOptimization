# NetworkOptimization
This project works on network optimization. Basically, we optimize data transmission, aiming particularly at Wireless Ad Hoc Networks (aka WANET or MANET). It's not yet established if this is only a research-project or if this code is going to be run on real network devices.

### The problem
As said above, the project optimizes data transmission in WANETs that are deployed in difficult terrain like mountains or marshes or any places affected by natural disasters. This type of networks is commonly used by search and rescue teams.
The very network on the ground is modeled as a graph on a lattice (the best lattice is hexagonal since GSM network uses it, but different approaches are also tested here), so the optimization quality is measured using different graph metrics, right now graph diameter is used - essentially, the less the graph diameter is, the faster data packets travel from the source node to the destination node in the network.
The network gets optimized by adding new nodes to it, which means that the job of the program is to figure out the best position for each of those new nodes.
A game-theoretic approach was proposed to solve this problem so later on in this description you will encounter words like "players", "strategies" and "payoffs".
The whole network consists of a number (usually a small number like 2 or 3) of subnets. Information can only be transmitted between the nodes that belong to the same subnet. So mathematically this can be seen as a graph consisting of a number of subgraphs whose configuration models the subnets. These subnets are supposed to be run by different network participants that are best though of as players. Each player has its own drone that has a data-transmitting device on it and can be deployed at any network position needed. Unlike other agents, drones can interact (sent data to and receive data from) with any node from any subnet which means a player can improve the network of another player. The strategy of each player is the position that they placed their drone at. The payoff of each player is calculated as the amount that the diameter of their subgraph gets reduced by after all the drones were deployed. And now we have a notion of a game here. This game is cooperative, which means that players do not tend to increase their own payoffs, but instead they are trying to maximize their collective payoff which is calculated as a sum of payoffs of each player (this might be a weighted sum or the collective payoff may be calculated using other approaches, but right now plain sum is used).

### Origins
Here's the background behind this work.
Initial problem was to figure out how to optimize data transmission in WANETs and write a program that finds out the solution for a given network configuration.
The first proposed solution was based on finding Nash equilibrium in a game with two participants (players) and it was discovered that Nash equilibrium does not always exist for this type of game. Only integer lattice (rectangular shape) was considered for simplicity.
Two following researches were also based on finding Nash equilibrium, but they tested two additional shaped of lattices - hexagonal and triangular.
The next research proposed maximizing collective payoff instead of finding out Nash equilibrium. The collective payoff was described above. This approach yields a set of solutions that are equally good at maximizing the collective payoff.
The most recent research done by [1](#dev-1) does the same stuff calculating collective payoff and figures out the set of solutions. What it does next is it defines a method to pick a single solution from the set so that it benefits the whole system in the best way possible (or better said the best way we came up with). The program in the repo is based on this research. It works on rectangular lattice and uses PageRank as a method for picking a solution.

### Researchers and Developers
<a name="dev-1">[1]</a> Sergey K., student at Saint-Petersburg State University (SPbSU), Saint-Petersburg, Russia. ssnnp@mail.ru, st054625@student.spbu.ru

## Project requirements
This project requires <strong>libpng</strong> and <strong>zlib</strong> to run properly. Those libraries are responsible for creating PNG images for visualization of the final results.