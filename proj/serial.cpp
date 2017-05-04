#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include <queue>
#include <utility>

#define MAX_LINE_WIDTH 100
#define INF 65535

using namespace std;
typedef pair<int, int> int_pair;

class AdjListNode {
public:
    AdjListNode(int node, int weight) : node(node), weight(weight) {}
    int node;
    int weight;
};

class AdjListEntry {
public:
    AdjListEntry(int source) : source(source) {}
    int source;
    vector<AdjListNode> neighbors;
};

void moore(int source, int nNodes, int dist[], vector<AdjListEntry>& adjList) {
    bool in_queue[nNodes + 1];
    for (int i = 1; i <= nNodes; ++ i) {
        dist[i] = INF;
        in_queue[i] = false;
    }
    dist[source] = 0;

    queue<int> fifo_q;
    fifo_q.push(source);
    in_queue[source] = true;

    int vi, vj, new_dist;
    while (fifo_q.size() > 0) {
        vi = fifo_q.front();
        in_queue[vi] = false;
        fifo_q.pop();

        for (int j = 0; j < adjList[vi].neighbors.size(); ++ j) {
            AdjListNode& adjNode = adjList[vi].neighbors[j];
            vj = adjNode.node;
            new_dist = adjNode.weight + dist[vi];

            if (new_dist < dist[vj] || dist[vj] == INF) {
                dist[vj] = new_dist;
                if (!in_queue[vj] && adjList[vj].neighbors.size() > 0) {
                    fifo_q.push(vj);
                    in_queue[vj] = true;
                }
            }
        }
    }
}

class VertexCompare {
    int *dist;
public:
    VertexCompare(int* dist) { this->dist = dist; }
    bool operator()(const int& v1, const int& v2) const {
        return dist[v1] > dist[v2];
    }
};

void dijkstra(int source, int nNodes, int dist[], vector<AdjListEntry>& adjList) {

    // priority_queue<int, vector<int>, VertexCompare> pq((VertexCompare(dist)));
    int visited[nNodes + 1];
    int unvisited = nNodes;

    for (int i = 1; i <= nNodes; ++ i) {
        visited[i] = false;
        dist[i] = (i == source) ? 0 : INF;
        // pq.push(i);
    }

    while (unvisited > 0) {
        // extract vertex u with minimum dist to source
        int u = 1, min_dist = INF;
        for (int i = 1; i <= nNodes; ++ i) {
            if (!visited[i] && dist[i] < min_dist) {
                min_dist = dist[i];
                u = i;
            }
        }
        visited[u] = true;
        unvisited -= 1;

        for (int i = 0; i < adjList[u].neighbors.size(); ++ i) {
            AdjListNode& adjNode = adjList[u].neighbors[i];
            int v = adjNode.node;
            if (visited[v]) continue;
            int new_dist = adjNode.weight + dist[u];
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
            }
        }
    }
}

int main(const int argc, const char** argv) {
    FILE *fp;
    if (argc <= 1 || (fp = fopen(argv[1], "r")) == NULL) {
        printf("\nArgument number error or file does not exist!\n");
        getchar();
        exit(1);
    }

    char strLine[MAX_LINE_WIDTH];
    char ch;
    int nNodes, nArcs;
    fgets(strLine, MAX_LINE_WIDTH, fp);
    fscanf(fp, "%c %s %d %d\n", &ch, strLine, &nNodes, &nArcs);

    vector<AdjListEntry> adjList;
    for (int i = 1; i <= nNodes; ++ i) {
        AdjListEntry newEntry(i);
        adjList.push_back(newEntry);
    }

    int src, dst, weight;
    for (int i = 0; i < nArcs; ++ i) {
        fscanf(fp, "%c %d %d %d\n", &ch, &src, &dst, &weight);
        AdjListNode newNode1(dst, weight);
        adjList[src].neighbors.push_back(newNode1);
        AdjListNode newNode2(src, weight);
        adjList[dst].neighbors.push_back(newNode2);
    }

    int dist[nNodes + 1];
    int source = 1;
    double tStart = omp_get_wtime(); // Start timing
    moore(source, nNodes, dist, adjList);
    double tEnd = omp_get_wtime(); // End timing

    for (int i = 1; i < nNodes; i += 1000) {
        printf("d(%d, %d) = %d\n", source, i, dist[i]);
    }
    printf("d(%d, %d) = %d\n", source, nNodes, dist[nNodes]);

//    tStart = omp_get_wtime(); // Start timing
//    dijkstra(source, nNodes, dist, adjList);
//    tEnd = omp_get_wtime(); // End timing
    printf("\nTook %.3lf s to compute.\n", tEnd - tStart);

    fclose(fp);
    return 0;
}
