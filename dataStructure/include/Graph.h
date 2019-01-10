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

    // Class to represent a graph
    class Graph {
        int V; // No. of vertices'

        // Pointer to an array containing adjacency listsList
        list<int> *adj;

        // A function used by topologicalSort
        void topologicalSortUtil(int v, bool visited[], stack<int> &Stack);

        std::set<int> nzB;

    public:

        Graph(int V, std::set<int> nzB); // Constructor

        // function to add an edge to graph
        void addEdge(int v, int w);

        // prints a Topological Sort of the complete graph
        void topologicalSort();

        void print();
    };

    Graph::Graph(int V, std::set<int> nzB) {
        this->nzB = nzB;
        this->V = V;
        adj = new list<int>[V];
    }

    void Graph::addEdge(int v, int w) {
        cout << v+1 << "->" << w+1 <<"     ";
        adj[v].push_back(w); // Add w to v’s list.
    }

// A recursive function used by topologicalSort
    void Graph::topologicalSortUtil(int v, bool visited[], stack<int> &Stack) {
        // Mark the current node as visited.
        visited[v] = true;

        // Recur for all the vertices adjacent to this vertex
        list<int>::iterator i;

        for (i = adj[v].begin(); i != adj[v].end(); ++i) {
            if (!visited[*i])
                topologicalSortUtil(*i, visited, Stack);
        }
        // Push current vertex to stack which stores result
        Stack.push(v);
    }

// The function to do Topological Sort. It uses recursive topologicalSortUtil()
    void Graph::topologicalSort() {

        stack<int> Stack;

        // Mark all the vertices as not visited
        bool *visited = new bool[V];
        for (int i = 0; i < V; i++)
            visited[i] = false;

        // Call the recursive helper function to store Topological Sort
        // starting from all vertices one by one
        for (int i:nzB)
            if (!visited[i])
                topologicalSortUtil(i, visited, Stack);

        cout << endl;
        // Print contents of stack
        while (!Stack.empty()) {
            cout << Stack.top() + 1 << " ";
            Stack.pop();
        }
    }
}

#endif //SPARSEOP2_GRAPH_H
