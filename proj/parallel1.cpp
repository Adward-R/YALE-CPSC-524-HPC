#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include <queue>
#include <list>

#define MAX_LINE_WIDTH 100
#define INF -1

using namespace std;

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

    while (fifo_q.size() > 0) {
        const int vi = fifo_q.front();
        in_queue[vi] = false;
        fifo_q.pop();

        #pragma omp parallel for
        for (int j = 0; j < adjList[vi].neighbors.size(); ++j) {
            const AdjListNode &adjNode = adjList[vi].neighbors[j];
            const int vj = adjNode.node;
            const int new_dist = adjNode.weight + dist[vi];

            if (new_dist < dist[vj] || dist[vj] == INF) {
                dist[vj] = new_dist;
                if (!in_queue[vj] && adjList[vj].neighbors.size() > 0) {
                    in_queue[vj] = true;
                    #pragma omp atomic capture
                    fifo_q.push(vj);
                }
            }
        }

    } // end of parallel section
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
    int nNodes;
    long nArcs;
    fgets(strLine, MAX_LINE_WIDTH, fp);
    fscanf(fp, "%c %s %d %ld\n", &ch, strLine, &nNodes, &nArcs);
    printf("nNodes = %d, nArcs = %ld\n", nNodes, nArcs);

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
    printf("\nTook %.3lf s to compute.\n", tEnd - tStart);

    fclose(fp);
    return 0;
}
