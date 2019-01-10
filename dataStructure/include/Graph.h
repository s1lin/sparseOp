//
// Created by shilei on 1/10/19.
//

#ifndef SPARSEOP2_GRAPH_H
#define SPARSEOP2_GRAPH_H

#include "Graph.h"
#include <iostream>
#include <cstdlib>
#include <list>
#include <stack>
#include <vector>

using namespace std;

namespace DataStructure {

    class Graph {

        // Number of vertices
        int numVer;

        //Number Of levels
        int level;

        // Pointer to an array containing adjacency lists
        list<int> *adj;

        // A function used by dfs
        void dfsUtil(int v, bool *visited, stack<int> &s);

        std::set<int> nzB;

    public:

        Graph(int numVer, std::set<int> nzB) {
            this->nzB = nzB;
            this->numVer = numVer;
            adj = new list<int>[numVer];
        }

        // function to add an edge to graph
        void addEdge(int v, int w) {
//            cout << v + 1 << "->" << w + 1 << "     ";
            adj[v].push_back(w); // Add w to v’s list.
        }

        // prints a dfs of the complete graph
        stack<int> dfs() {

            stack<int> s;
            level = 0;

            // Mark all the vertices as not visited
            bool *visited = new bool[numVer];
            for (int i = 0; i < numVer; i++)
                visited[i] = false;

            // Call the recursive helper function to store Topological Sort
            // starting from all vertices one by one
            for (int i:nzB)
                if (!visited[i])
                    dfsUtil(i, visited, s);

//            stack<int> s2 = s;
////             Print contents of stack
//            while (!s2.empty()) {
//                cout << s2.top() + 1 << " ";
//                s2.pop();
//            }

            return s;
        }

    };

    void Graph::dfsUtil(int v, bool *visited, stack<int> &s) {
        // Mark the current node as visited.
        visited[v] = true;

        // Recur for all the vertices adjacent to this vertex
        list<int>::iterator i;

        for (i = adj[v].begin(); i != adj[v].end(); ++i) {
            if (!visited[*i])
                dfsUtil(*i, visited, s);
        }
        // Push current vertex to stack which stores result
        s.push(v);

    }

}

#endif //SPARSEOP2_GRAPH_H
