/******************************************************************************
 * File: Trailblazer.cpp
 *
 * Implementation of the graph algorithms that comprise the Trailblazer
 * assignment.
 */

#include "Trailblazer.h"
#include "TrailblazerGraphics.h"
#include "TrailblazerTypes.h"
#include "TrailblazerPQueue.h"
#include "random.h"
#include "foreach.h"
#include <limits>
using namespace std;

/* Function: shortestPath
 *
 * Finds the shortest path between the locations given by start and end in the
 * specified world.	 The cost of moving from one edge to the next is specified
 * by the given cost function.	The resulting path is then returned as a
 * Vector<Loc> containing the locations to visit in the order in which they
 * would be visited.	If no path is found, this function should report an
 * error.
 *
 * In Part Two of this assignment, you will need to add an additional parameter
 * to this function that represents the heuristic to use while performing the
 * search.  Make sure to update both this implementation prototype and the
 * function prototype in Trailblazer.h.
 */
Vector<Loc>
shortestPath(Loc start,
    Loc end,
    Grid<double>& world,
    double costFn(Loc from, Loc to, Grid<double>& world),
    double heuristic(Loc start, Loc end, Grid<double>& world)) {
    // TODO: Fill this in!
    
    Grid<BookKeep> bk(world.numRows(), world.numCols());


    for (int i = 0; i < bk.numRows(); i++) {
        for (int j = 0; j < bk.numCols(); j++) {
            bk[i][j].color = GRAY;
            bk[i][j].parent = NIL;
        }
    }

    bk[start.row][start.col].color = YELLOW;
    bk[start.row][start.col].distance = heuristic(start, end, world);

    TrailblazerPQueue<Loc> queue;
    queue.enqueue(start, 0);

    Loc loc = start;

    while (true) {

        loc = queue.dequeueMin();

        bk[loc.row][loc.col].color = GREEN;
        colorCell(world, loc, GREEN);

        if (loc == end) break;

        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {

                Loc v = makeLoc(loc.row + i, loc.col + j);

                if (world.inBounds(v.row, v.col) && v != loc) {
                    if (bk[v.row][v.col].color == GRAY) {

                        bk[v.row][v.col].color = YELLOW;
                        colorCell(world, v, YELLOW);

                        double L = costFn(loc, v, world);
                        double dist = bk[v.row][v.col].distance;
                        bk[v.row][v.col].parent = loc;
                        bk[v.row][v.col].distance = dist + L;

                        queue.enqueue(v, dist + L + heuristic(v, end, world));
                    }

                    if (bk[v.row][v.col].color == YELLOW) {
                        double L = costFn(loc, v, world);
                        double dist = bk[loc.row][loc.col].distance;
                        if (bk[v.row][v.col].distance > L + dist) {

                            bk[v.row][v.col].distance = L + dist;
                            bk[v.row][v.col].parent = loc;

                            queue.decreaseKey(v, L + dist + heuristic(v, end, world));
                        }
                    }
                }
            }
        }
    }

    Vector<Loc> reversePath;

    while (end != NIL) {
        reversePath.push_back(end);
        end = bk[end.row][end.col].parent;
    }

    Vector<Loc> path;
    for (int i = reversePath.size() - 1; i >= 0; i--) {
        path.push_back(reversePath[i]);
    }

    return path;
}

struct Cluster {
    int parent;
    int rank;
};

int find(int node, Vector<Cluster>& clusters) {
    if (clusters[node].parent != node) {
        clusters[node].parent = find(clusters[node].parent, clusters);
    }
    return clusters[node].parent;
}

void shuffle(Vector<Edge>& vec) {
    int n = vec.size();
    for (int i = n - 1; i > 0; i--) {

        int j = randomInteger(0, i);


        Edge temp = vec[i];
        vec[i] = vec[j];
        vec[j] = temp;
    }
}


void unionNodes(int node1, int node2, Vector<Cluster>& clusters) {
    int rootA = find(node1, clusters);
    int rootB = find(node2, clusters);

    if (rootB != rootA) {
        if (clusters[rootA].rank > clusters[rootA].rank) {
            clusters[rootB].parent = rootA;
        }
        else if (clusters[rootA].rank < clusters[rootB].rank) {

            clusters[rootA].parent = rootB;

        }
        else {
            clusters[rootB].parent = rootA;
            clusters[rootA].rank++;
        }
    }
}

Set<Edge> createMaze(int numRows, int numCols) {
    int sq = numRows * numCols;

    Set<Edge> mazeEdges;
    Vector<Cluster> clusters(sq);

  
    for (int i = 0; i < sq; ++i) {
        clusters[i].parent = i;
        clusters[i].rank = 0;
    }


    Vector<Edge> allEdges;
    for (int row = 0; row < numRows; ++row) {
        for (int col = 0; col < numCols; ++col) {
            Loc current = makeLoc(row, col);

            //make horizontal
            if (col + 1 < numCols) {
                Loc right = makeLoc(row, col + 1);
                allEdges.add(makeEdge(current, right));
            }

            //now vertical 
            if (row + 1 < numRows) {
                Loc down = makeLoc(row + 1, col);
                allEdges.add(makeEdge(current, down));
            }
        }
    }

    shuffle(allEdges);

    for (Edge edges : allEdges) {
        int startNode = edges.start.row * numCols + edges.start.col;
        int endNode = edges.end.row * numCols + edges.end.col;

        if (find(startNode, clusters) != find(endNode, clusters)) {

            unionNodes(startNode, endNode, clusters);
            mazeEdges.add(edges);
        }

        if (mazeEdges.size() == (sq - 1)) break;
    }
    return mazeEdges;
}
