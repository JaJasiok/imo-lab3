#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <climits>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <numeric>
#include <variant>
#include <chrono>

using namespace std;

vector<vector<int>> readKroaFile(const string &filename)
{
    ifstream file(filename);
    string line;
    vector<vector<int>> verticesCoords;

    if (file.is_open())
    {
        // Skip the header lines
        for (int i = 0; i < 7; i++)
        {
            getline(file, line);
        }

        // Read the coordinates and populate the cost matrix
        do
        {
            istringstream iss(line);
            int index, x, y;
            iss >> index >> x >> y;
            verticesCoords.push_back({x, y});
            // cout << index << " " << x << " " << y << endl;
        } while (getline(file, line) && line != "EOF");

        // verticesCoords.pop_back();

        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }

    return verticesCoords;
}

vector<vector<int>> createDistanceMatrix(const vector<vector<int>> &verticesCoords)
{
    int numVertices = verticesCoords.size();
    vector<vector<int>> distanceMatrix(numVertices, vector<int>(numVertices));

    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            double x1 = verticesCoords[i][0];
            double y1 = verticesCoords[i][1];
            double x2 = verticesCoords[j][0];
            double y2 = verticesCoords[j][1];

            double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
            distanceMatrix[i][j] = int(round(distance));

            // cout << i << " " << j << " " << distanceMatrix[i][j] << endl;
        }
    }

    return distanceMatrix;
}

class Vertex
{
public:
    int id;
    Vertex *prev;
    Vertex *next;

    Vertex(int _id) : id(_id), prev(nullptr), next(nullptr) {}
};

class Edge
{
public:
    Vertex *src;
    Vertex *dest;
    int distance;

    Edge(Vertex *_src, Vertex *_dest, int _distance) : src(_src), dest(_dest), distance(_distance) {}

    void remove()
    {
        delete this;
    }
};

class Graph
{
public:
    vector<Vertex *> vertices;
    vector<Edge *> edges;
    int distance = 0;

    void addVertex(Vertex *v)
    {
        if (v)
        {
            if (find(vertices.begin(), vertices.end(), v) == vertices.end())
            {
                vertices.push_back(v);
            }
        }
    }

    void addEdge(Vertex *src, Vertex *dest, int distance = 0)
    {
        try
        {
            if (src && dest && src != dest)
            {
                for (Edge *e : edges)
                {
                    if ((e->src == src && e->dest == dest) || (e->src == dest && e->dest == src))
                    {
                        return;
                    }
                }

                Edge *e = new Edge(src, dest, distance);
                src->next = dest;
                dest->prev = src;
                this->distance += distance;
                edges.push_back(e);
            }
        }
        catch (exception &e)
        {
            cout << e.what() << endl;
        }
    }

    void removeVertex(Vertex *v)
    {
        if (v)
        {
            removeEdge(v->prev, v);
            removeEdge(v, v->next);

            vertices.erase(remove(vertices.begin(), vertices.end(), v), vertices.end());
        }
    }

    void removeEdge(Vertex *src, Vertex *dest)
    {
        if (src && dest)
        {
            for (Edge *e : edges)
            {
                if (e->src == src && e->dest == dest)
                {
                    this->distance -= e->distance;
                    edges.erase(remove(edges.begin(), edges.end(), e), edges.end());
                    // e->remove();
                    break;
                }
            }
            src->next = nullptr;
            dest->prev = nullptr;
        }
    }

    Vertex *findVertex(int id)
    {
        for (Vertex *v : vertices)
        {
            if (v->id == id)
                return v;
        }
        return nullptr;
    }

    Edge *findEdge(int src_id, int dest_id)
    {
        for (Edge *e : edges)
        {
            if (e->src->id == src_id && e->dest->id == dest_id)
            {
                return e;
            }
        }
        return nullptr;
    }

    vector<Edge *> findPath(Vertex *start, Vertex *end)
    {
        vector<Edge *> path;
        Vertex *currentVertex = start;

        do
        {
            for (Edge *e : edges)
            {
                if (e->src == currentVertex)
                {
                    path.push_back(e);
                    currentVertex = e->dest;
                    break;
                }
            }
        } while (currentVertex != end);

        return path;
    }
};

class Move
{
public:
    pair<Vertex *, Vertex *> vertices;
    pair<Edge *, Edge *> edges;
    Graph *graph;
    int delta;

    Move(pair<Vertex *, Vertex *> _vertices = {nullptr, nullptr}, pair<Edge *, Edge *> _edges = {nullptr, nullptr}, Graph *_graph = nullptr, int _delta = 0)
        : vertices(_vertices), edges(_edges), graph(_graph), delta(_delta) {}
};

void saveGraphs(const vector<Graph> &graphs, const string &filename)
{
    ofstream file(filename);
    if (file.is_open())
    {
        for (Graph g : graphs)
        {
            for (Edge *e : g.edges)
            {
                file << e->src->id << " " << e->dest->id << endl;
            }
            file << endl;
        }
        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }
}

vector<Graph> randomCycles(const vector<vector<int>> &distanceMatrix)
{
    int numVertices = distanceMatrix.size();

    std::vector<int> numbers(numVertices);
    iota(numbers.begin(), numbers.end(), 0);

    random_shuffle(numbers.begin(), numbers.end());

    vector<Graph> cycles(2);

    std::vector<int> group1(numbers.begin(), numbers.begin() + numVertices / 2);
    std::vector<int> group2(numbers.begin() + numVertices / 2, numbers.end());

    for (int i : group1)
    {
        Vertex *v = new Vertex(i);
        cycles[0].addVertex(v);
    }

    for (int i : group2)
    {
        Vertex *v = new Vertex(i);
        cycles[1].addVertex(v);
    }

    for (int i = 0; i < numVertices / 2; i++)
    {
        if (i < (numVertices / 2) - 1)
        {
            Vertex *v1 = cycles[0].vertices[i];
            Vertex *v2 = cycles[0].vertices[i + 1];
            // cout << v1->id << " " << v2->id << " " << distanceMatrix[v1->id][v2->id] << endl;
            cycles[0].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
        else
        {
            Vertex *v1 = cycles[0].vertices[i];
            Vertex *v2 = cycles[0].vertices[0];
            // cout << v1->id << " " << v2->id << " " << distanceMatrix[v1->id][v2->id] << endl;
            cycles[0].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
    }

    for (int i = 0; i < numVertices / 2; i++)
    {
        if (i < (numVertices / 2) - 1)
        {
            Vertex *v1 = cycles[1].vertices[i];
            Vertex *v2 = cycles[1].vertices[i + 1];
            // cout << distanceMatrix[v1->id][v2->id] << endl;
            cycles[1].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
        else
        {
            Vertex *v1 = cycles[1].vertices[i];
            Vertex *v2 = cycles[1].vertices[0];
            // cout << distanceMatrix[v1->id][v2->id] << endl;
            cycles[1].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
    }

    return cycles;
}

vector<Graph> greedyCycles(const vector<vector<int>> &distanceMatrix, int startId)
{
    int numVertices = distanceMatrix.size();
    vector<Graph> cycles(2);
    vector<bool> visited(numVertices, false);

    // Choose the first vertex randomly
    Vertex *startVertex1 = new Vertex(startId);
    cycles[0].addVertex(startVertex1);
    visited[startVertex1->id] = true;

    // Find the furthest vertex from startVertex1
    int furthestVertex = -1;
    int maxDistance = -1;
    for (int j = 0; j < numVertices; j++)
    {
        if (!visited[j] && distanceMatrix[startVertex1->id][j] > maxDistance)
        {
            maxDistance = distanceMatrix[startVertex1->id][j];
            furthestVertex = j;
        }
    }

    Vertex *startVertex2 = new Vertex(furthestVertex);
    cycles[1].addVertex(startVertex2);
    visited[startVertex2->id] = true;

    // Find the nearest neighbour for each starting vertex
    for (Graph &cycle : cycles)
    {
        Vertex *vertex = cycle.vertices.front();
        int vertexId = vertex->id;
        int minDistance = INT_MAX;
        int nearestNeighbourId = -1;

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                if (distanceMatrix[vertexId][j] < minDistance)
                {
                    minDistance = distanceMatrix[vertexId][j];
                    nearestNeighbourId = j;
                }
            }
        }

        Vertex *nearestNeighbour = new Vertex(nearestNeighbourId);

        cycle.addVertex(nearestNeighbour);
        cycle.addEdge(vertex, nearestNeighbour, minDistance);
        cycle.addEdge(nearestNeighbour, vertex);
        visited[nearestNeighbourId] = true;
    }

    // Building the rest of the cycle
    for (int i = 0; i < numVertices - 4; i++)
    {
        int minDistance = INT_MAX;
        int vertexId = -1;
        pair<Edge *, Graph *> minPair;

        vector<pair<Edge *, Graph *>> edgesInGraphs;
        if (cycles[0].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[1].edges)
            {
                edgesInGraphs.push_back({e, &cycles[1]});
            }
        }
        else if (cycles[1].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[0].edges)
            {
                edgesInGraphs.push_back({e, &cycles[0]});
            }
        }
        else
        {
            for (Graph &cycle : cycles)
            {
                for (Edge *e : cycle.edges)
                {
                    edgesInGraphs.push_back({e, &cycle});
                }
            }
        }

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                for (pair<Edge *, Graph *> singlePair : edgesInGraphs)
                {
                    Edge *e = singlePair.first;
                    int distanceSum = cycles[0].distance + cycles[1].distance;
                    if (distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j] < minDistance)
                    {
                        minDistance = distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j];
                        vertexId = j;
                        minPair = singlePair;
                    }
                }
            }
        }

        Graph *cyclePtr = minPair.second;
        Graph &cycle = *cyclePtr;

        Edge *minEdge = minPair.first;

        Vertex *newVertex = new Vertex(vertexId);
        cycle.addVertex(newVertex);

        // if (cycle.vertices.size() > 3)
        // {
        cycle.removeEdge(minEdge->src, minEdge->dest);
        // }

        cycle.addEdge(minEdge->src, newVertex, distanceMatrix[minEdge->src->id][vertexId]);
        cycle.addEdge(newVertex, minEdge->dest, distanceMatrix[newVertex->id][minEdge->dest->id]);

        if (cycle.vertices.size() == 3)
        {
            cycle.addEdge(minEdge->dest, minEdge->src, minEdge->distance);
        }

        visited[vertexId] = true;
    }

    return cycles;
}

void swapVerticesBetweenCycles(Vertex *vertex1, Vertex *vertex2, vector<Graph> &graphs, const vector<vector<int>> &distanceMatrix)
{
    Graph *graph1;
    Graph *graph2;

    if(graphs[0].findVertex(vertex1->id) != nullptr)
    {
        graph1 = &graphs[0];
        graph2 = &graphs[1];
    }
    else
    {
        graph1 = &graphs[1];
        graph2 = &graphs[0];
    }


    Vertex *vertex1Prev = vertex1->prev;
    Vertex *vertex1Next = vertex1->next;

    Vertex *vertex2Prev = vertex2->prev;
    Vertex *vertex2Next = vertex2->next;

    graph1->removeVertex(vertex1);

    graph2->removeVertex(vertex2);

    graph1->addVertex(vertex2);
    graph1->addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
    graph1->addEdge(vertex2, vertex1Next, distanceMatrix[vertex2->id][vertex1Next->id]);

    graph2->addVertex(vertex1);
    graph2->addEdge(vertex2Prev, vertex1, distanceMatrix[vertex2Prev->id][vertex1->id]);
    graph2->addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
}

void swapEdgesInGraph(Edge *edge1, Edge *edge2, Graph *graph, const vector<vector<int>> &distanceMatrix)
{
    Vertex *vertex11 = edge1->src;
    Vertex *vertex12 = edge1->dest;
    Vertex *vertex21 = edge2->src;
    Vertex *vertex22 = edge2->dest;

    // for(Edge *e: graph->edges)
    // {
    //     cout << e->src->id << " " << e->dest->id << endl;
    // }

    graph->removeEdge(vertex11, vertex12);
    graph->removeEdge(vertex21, vertex22);

    vector<Edge *> pathToReverse = graph->findPath(vertex12, vertex21);

    for (Edge *edge : pathToReverse)
    {
        swap(edge->src, edge->dest);
        swap(edge->dest->next, edge->dest->prev);
        if (edge->src->next == nullptr)
        {
            swap(edge->src->next, edge->src->prev);
        }
    }

    graph->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
    graph->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
}

void swapVerticesInCycle(Vertex *vertex1, Vertex *vertex2, Graph *graph, const vector<vector<int>> &distanceMatrix)
{
    Vertex *vertex1Prev = vertex1->prev;
    Vertex *vertex1Next = vertex1->next;

    Vertex *vertex2Prev = vertex2->prev;
    Vertex *vertex2Next = vertex2->next;


    graph->removeEdge(vertex1Prev, vertex1);
    graph->removeEdge(vertex1, vertex1Next);
    graph->removeEdge(vertex2Prev, vertex2);
    graph->removeEdge(vertex2, vertex2Next);

    graph->addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
    graph->addEdge(vertex2, vertex1Next, distanceMatrix[vertex2->id][vertex1Next->id]);
    graph->addEdge(vertex2Prev, vertex1, distanceMatrix[vertex2Prev->id][vertex1->id]);
    graph->addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
}

vector<Graph> steepestLocalSearch(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    int bestDelta = 0;

    int i = 0;

    do
    {
        vector<Move *> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                Move *move = new Move({vertex1, vertex2}, {nullptr, nullptr}, nullptr, 0);
                moves.push_back(move);
            }
        }

        for (Graph &graph : cycles)
        {
            for (Edge *edge1 : graph.edges)
            {
                for (Edge *edge2 : graph.edges)
                {
                    if (edge1 != edge2)
                    {
                        Vertex *vertex11 = edge1->src;
                        Vertex *vertex12 = edge1->dest;
                        Vertex *vertex21 = edge2->src;
                        Vertex *vertex22 = edge2->dest;

                        if (vertex11->id != vertex12->id && vertex11->id != vertex21->id && vertex11->id != vertex22->id && vertex12->id != vertex21->id && vertex12->id != vertex22->id && vertex21->id != vertex22->id)
                        {
                            Move *move = new Move({nullptr, nullptr}, {edge1, edge2}, &graph, 0);
                            moves.push_back(move);
                        }
                    }
                }
            }
        }

        Move *bestMove = nullptr;

        bestDelta = 0;

        for (auto move : moves)
        {
            int delta = 0;

            if (move->graph == nullptr)
            {
                Vertex *vertex1 = move->vertices.first;
                Vertex *vertex2 = move->vertices.second;

                delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];
            }
            else
            {
                Edge *edge1 = move->edges.first;
                Edge *edge2 = move->edges.second;
                Graph *graphToSwap = move->graph;

                Vertex *vertex11 = edge1->src;
                Vertex *vertex12 = edge1->dest;
                Vertex *vertex21 = edge2->src;
                Vertex *vertex22 = edge2->dest;

                delta = distanceMatrix[vertex11->id][vertex21->id] + distanceMatrix[vertex12->id][vertex22->id] - edge1->distance - edge2->distance;
            }

            if (delta < bestDelta)
            {
                bestMove = move;
                bestDelta = delta;
            }
        }

        if (bestDelta < 0)
        {
            if (bestMove->graph == nullptr)
            {
                Vertex *vertex1 = bestMove->vertices.first;
                Vertex *vertex2 = bestMove->vertices.second;

                swapVerticesBetweenCycles(vertex1, vertex2, cycles, distanceMatrix);
            }
            else
            {
                Edge *edge1 = bestMove->edges.first;
                Edge *edge2 = bestMove->edges.second;

                Graph *graphToSwap = bestMove->graph;

                swapEdgesInGraph(edge1, edge2, graphToSwap, distanceMatrix);
            }
        }

        i++;
    } while (bestDelta < 0);

    return cycles;
}

vector<Graph> localSearchWithHistory(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    // int it = 0;

    int i = 0;

    vector<Move *> moves;

    for (Vertex *vertex1 : cycles[0].vertices)
    {
        for (Vertex *vertex2 : cycles[1].vertices)
        {
            int delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];

            if (delta < 0)
            {
                Move *move = new Move({vertex1, vertex2}, {nullptr, nullptr}, nullptr, delta);
                moves.push_back(move);
            }
        }
    }

    for (Graph &graph : cycles)
    {
        for (Edge *edge1 : graph.edges)
        {
            for (Edge *edge2 : graph.edges)
            {
                if (edge1 != edge2)
                {
                    Vertex *vertex11 = edge1->src;
                    Vertex *vertex12 = edge1->dest;
                    Vertex *vertex21 = edge2->src;
                    Vertex *vertex22 = edge2->dest;

                    if (vertex11->id != vertex12->id && vertex11->id != vertex21->id && vertex11->id != vertex22->id && vertex12->id != vertex21->id && vertex12->id != vertex22->id && vertex21->id != vertex22->id)
                    {
                        int delta = distanceMatrix[vertex11->id][vertex21->id] + distanceMatrix[vertex12->id][vertex22->id] - edge1->distance - edge2->distance;
                        if (delta < 0)
                        {
                            Move *move = new Move({nullptr, nullptr}, {edge1, edge2}, &graph, delta);
                            moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    sort(moves.begin(), moves.end(), [](Move *a, Move *b)
         { return a->delta < b->delta; });

    do
    {
        for (i = 0; i < moves.size(); i++)
        {
            auto move = moves[i];

            if (move->graph == nullptr)
            {
                Vertex *vertex1 = move->vertices.first;
                Vertex *vertex2 = move->vertices.second;

                int delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];

                move->delta = delta;
            }
            else
            {
                Edge *edge1 = move->edges.first;
                Edge *edge2 = move->edges.second;
                Graph *graphToSwap = move->graph;

                Vertex *vertex11 = edge1->src;
                Vertex *vertex12 = edge1->dest;
                Vertex *vertex21 = edge2->src;
                Vertex *vertex22 = edge2->dest;

                int delta = distanceMatrix[vertex11->id][vertex21->id] + distanceMatrix[vertex12->id][vertex22->id] - edge1->distance - edge2->distance;

                move->delta = delta;
            }
        }

        sort(moves.begin(), moves.end(), [](Move *a, Move *b)
             { return a->delta < b->delta; });

        for (i = 0; i < moves.size(); i++)
        {
            // cout << it << " " << i << "/" << moves.size() << endl;

            auto move = moves[i];

            if (move->graph == nullptr)
            {
                Vertex *vertex1 = move->vertices.first;
                Vertex *vertex2 = move->vertices.second;

                if (cycles[0].findVertex(vertex1->id) == nullptr || cycles[1].findVertex(vertex2->id) == nullptr)
                {
                    moves.erase(moves.begin() + i);
                    i--;
                    continue;
                }

                swapVerticesBetweenCycles(vertex1, vertex2, cycles, distanceMatrix);

                for (Vertex *vertex : cycles[1].vertices)
                {
                    int delta = distanceMatrix[vertex->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex->next->id] + distanceMatrix[vertex2->prev->id][vertex->id] + distanceMatrix[vertex->id][vertex2->next->id] - distanceMatrix[vertex->prev->id][vertex->id] - distanceMatrix[vertex->id][vertex->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];

                    if (delta < 0)
                    {
                        Move *move = new Move({vertex2, vertex}, {nullptr, nullptr}, nullptr, delta);
                        moves.push_back(move);
                    }
                }
                for (Vertex *vertex : cycles[0].vertices)
                {
                    int delta = distanceMatrix[vertex->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex->next->id] + distanceMatrix[vertex1->prev->id][vertex->id] + distanceMatrix[vertex->id][vertex1->next->id] - distanceMatrix[vertex->prev->id][vertex->id] - distanceMatrix[vertex->id][vertex->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id];

                    if (delta < 0)
                    {
                        Move *move = new Move({vertex, vertex1}, {nullptr, nullptr}, nullptr, delta);
                        moves.push_back(move);
                    }
                }

                // it++;
                break;
            }
            else
            {
                Edge *edge1 = move->edges.first;
                Edge *edge2 = move->edges.second;
                Graph *graphToSwap = move->graph;

                Vertex *vertex11 = edge1->src;
                Vertex *vertex12 = edge1->dest;
                Vertex *vertex21 = edge2->src;
                Vertex *vertex22 = edge2->dest;

                if ((graphToSwap->findEdge(vertex11->id, vertex12->id) != nullptr && graphToSwap->findEdge(vertex21->id, vertex22->id) != nullptr))
                {
                    swapEdgesInGraph(edge1, edge2, graphToSwap, distanceMatrix);

                    for (Edge *edge : graphToSwap->edges)
                    {
                        if (edge != edge1)
                        {
                            Vertex *vertex1 = edge->src;
                            Vertex *vertex2 = edge->dest;

                            if (vertex11->id != vertex12->id && vertex11->id != vertex1->id && vertex11->id != vertex2->id && vertex12->id != vertex1->id && vertex12->id != vertex2->id && vertex1->id != vertex2->id)
                            {
                                int delta = distanceMatrix[edge->src->id][vertex11->id] + distanceMatrix[edge->dest->id][vertex12->id] - edge->distance - edge1->distance;
                                if (delta < 0)
                                {
                                    Move *move1 = new Move({nullptr, nullptr}, {edge, edge1}, graphToSwap, delta);
                                    Move *move2 = new Move({nullptr, nullptr}, {edge1, edge}, graphToSwap, delta);

                                    moves.push_back(move1);
                                    moves.push_back(move2);
                                }
                            }
                        }
                        if (edge != edge2)
                        {
                            Vertex *vertex1 = edge->src;
                            Vertex *vertex2 = edge->dest;

                            if (vertex1->id != vertex2->id && vertex1->id != vertex21->id && vertex1->id != vertex22->id && vertex2->id != vertex21->id && vertex2->id != vertex22->id && vertex21->id != vertex22->id)
                            {
                                int delta = distanceMatrix[edge->src->id][vertex21->id] + distanceMatrix[edge->dest->id][vertex22->id] - edge->distance - edge2->distance;
                                if (delta < 0)
                                {
                                    Move *move1 = new Move({nullptr, nullptr}, {edge, edge2}, graphToSwap, delta);
                                    Move *move2 = new Move({nullptr, nullptr}, {edge2, edge}, graphToSwap, delta);

                                    moves.push_back(move1);
                                    moves.push_back(move2);
                                }
                            }
                        }
                    }

                    // it++;
                    break;
                }

                if ((graphToSwap->findEdge(vertex11->id, vertex12->id) != nullptr && graphToSwap->findEdge(vertex22->id, vertex21->id) != nullptr) || (graphToSwap->findEdge(vertex12->id, vertex11->id) != nullptr && graphToSwap->findEdge(vertex21->id, vertex22->id) != nullptr))
                {
                    continue;
                }


                moves.erase(moves.begin() + i);
                i--;
                continue;
            }
        }

    } while (moves[i]->delta < 0);

    return cycles;
}

vector<Graph> localSearchWithCandidates(vector<Graph> &cycles, const vector<vector<int>> &distanceMatrix)
{
    vector<Move *> moves;

    do
    {
        moves.clear();

        for (int i = 0; i < cycles.size(); i++)
        {
            Graph &graph = cycles[i];

            for (Vertex *vertex1 : graph.vertices)
            {
                vector<pair<int, int>> closestVertices;
                for (int j = 0; j < vertex1->id; j++)
                {
                    if (vertex1->id != j)
                    {
                        int vertex2id = j;
                        closestVertices.push_back({distanceMatrix[vertex1->id][j], vertex2id});
                    }
                }

                sort(closestVertices.begin(), closestVertices.end());
                if(closestVertices.size() > 20)
                {
                    closestVertices.resize(20);
                }

                for (auto closestVertex : closestVertices)
                {

                    if (graph.findVertex(closestVertex.second) == nullptr)
                    {
                        Vertex* vertex2;
                        if(i == 0){
                            vertex2 = cycles[1].findVertex(closestVertex.second);
                        }
                        else{
                            vertex2 = cycles[0].findVertex(closestVertex.second);
                        }

                        int delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];

                        if (delta < 0)
                        {
                            Move *move = new Move({vertex1, vertex2}, {nullptr, nullptr}, nullptr, delta);
                            moves.push_back(move);
                        }
                    }
                    else
                    {
                        Vertex *vertex2 = graph.findVertex(closestVertex.second);
                        if (graph.findEdge(vertex1->id, vertex2->id) != nullptr && graph.findEdge(vertex2->id, vertex1->id) != nullptr)
                        {
                            int delta = 0;

                            Edge* edge1 = graph.findEdge(vertex1->prev->id, vertex1->id);
                            Edge* edge2 = graph.findEdge(vertex2->prev->id, vertex2->id);

                            delta = distanceMatrix[vertex1->prev->id][vertex2->prev->id] + distanceMatrix[vertex1->id][vertex2->id] - edge1->distance - edge2->distance;

                            if (delta < 0)
                            {
                                Move *move = new Move({nullptr, nullptr}, {edge1, edge2}, &graph, delta);
                                moves.push_back(move);
                            }

                            Edge* edge3 = graph.findEdge(vertex1->id, vertex1->next->id);
                            Edge* edge4 = graph.findEdge(vertex2->id, vertex2->next->id);

                            delta = distanceMatrix[vertex1->id][vertex2->id] + distanceMatrix[vertex1->next->id][vertex2->next->id] - edge3->distance - edge4->distance;

                            if (delta < 0)
                            {
                                Move *move = new Move({nullptr, nullptr}, {edge3, edge4}, &graph, delta);
                                moves.push_back(move);
                            }
                        }
                    }
                }
            }
        }

        sort(moves.begin(), moves.end(), [](Move *a, Move *b)
             { return a->delta < b->delta; });

        Move *bestMove = moves[0];

        if (bestMove->graph == nullptr)
        {
            Vertex *vertex1 = bestMove->vertices.first;
            Vertex *vertex2 = bestMove->vertices.second;

            swapVerticesBetweenCycles(vertex1, vertex2, cycles, distanceMatrix);
        }
        else
        {
            Edge *edge1 = bestMove->edges.first;
            Edge *edge2 = bestMove->edges.second;

            Graph *graphToSwap = bestMove->graph;

            swapEdgesInGraph(edge1, edge2, graphToSwap, distanceMatrix);
        }

    } while (moves.size() > 0);

    return cycles;
}

int main()
{
    // srand(time(0));

    vector<vector<int>> verticesCoords = readKroaFile("kroB200.tsp");

    vector<vector<int>> distanceMatrix = createDistanceMatrix(verticesCoords);

    pair<int, int> bestValue = {INT_MAX, -1};
    pair<int, int> worstValue = {-1, -1};
    long averageValue = 0;

    pair<int, int> bestTime = {INT_MAX, -1};
    pair<int, int> worstTime = {-1, -1};
    long averageTime = 0;

    vector<Graph> bestResult;

    int n = 100;

    for (int i = 0; i < n; i++)
    {
        vector<Graph> cycles = randomCycles(distanceMatrix);

        chrono::steady_clock::time_point begin = chrono::steady_clock::now();

        vector<Graph> result = localSearchWithHistory(cycles, distanceMatrix);

        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        int elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();

        cout << result[0].distance + result[1].distance << endl;

        if (result[0].distance + result[1].distance < bestValue.first)
        {
            bestValue = {result[0].distance + result[1].distance, i};
            bestResult = result;
        }

        if (result[0].distance + result[1].distance > worstValue.first)
        {
            worstValue = {result[0].distance + result[1].distance, i};
        }

        averageValue += result[0].distance + result[1].distance;

        if (elapsed > worstTime.first)
        {
            worstTime = {elapsed, i};
        }

        if (elapsed < bestTime.first)
        {
            bestTime = {elapsed, i};
        }

        averageTime += elapsed;

        cout << i << endl;
    }

    averageValue /= n;
    averageTime /= n;

    cout << "localSearchWithHistory" << endl;
    cout << averageValue << " (" << bestValue.first << " – " << worstValue.first << ")" << endl;
    cout << averageTime << " (" << bestTime.first << " – " << worstTime.first << ")" << endl;

    // saveGraphs(bestResult, "steepestLocalSearchB.txt");

    return 0;
}