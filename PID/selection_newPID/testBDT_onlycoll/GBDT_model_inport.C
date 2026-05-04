#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>

struct Node { int left, right, feature; double threshold, value; };

struct Tree {
    std::vector<Node> nodes;
    double predict(const std::vector<double>& x) const {
        int cur = 0;
        while (nodes[cur].feature != -2)
            cur = x[nodes[cur].feature] <= nodes[cur].threshold
                  ? nodes[cur].left : nodes[cur].right;
        return nodes[cur].value;
    }
};

struct GBDTModel {
    double learning_rate;
    int n_classes;
    std::vector<double> init_scores;  // uno per classe
    std::vector<Tree> trees;          // n_stages * n_classes alberi

    // Ritorna il vettore di probabilità per ogni classe (softmax)
    std::vector<double> predict_proba(const std::vector<double>& x) const {
        int n_stages = trees.size() / n_classes;

        // Accumula score per ogni classe
        std::vector<double> scores = init_scores;
        for (int s = 0; s < n_stages; s++)
            for (int k = 0; k < n_classes; k++)
                scores[k] += learning_rate * trees[s * n_classes + k].predict(x);

        // Softmax
        double max_s = *std::max_element(scores.begin(), scores.end());
        std::vector<double> proba(n_classes);
        double sum = 0.0;
        for (int k = 0; k < n_classes; k++) {
            proba[k] = std::exp(scores[k] - max_s);
            sum += proba[k];
        }
        for (auto& p : proba) p /= sum;
        return proba;
    }

    int predict_class(const std::vector<double>& x) const {
        auto proba = predict_proba(x);
        return std::max_element(proba.begin(), proba.end()) - proba.begin();
    }
};

GBDTModel load_model(const char* filename) {
    std::ifstream f(filename);
    GBDTModel model;

    int n_stages;
    f >> model.learning_rate >> n_stages >> model.n_classes;

    model.init_scores.resize(model.n_classes);
    for (int k = 0; k < model.n_classes; k++)
        f >> model.init_scores[k];

    for (int s = 0; s < n_stages * model.n_classes; s++) {
        int n_nodes; f >> n_nodes;
        Tree tree;
        for (int i = 0; i < n_nodes; i++) {
            Node n;
            f >> n.left >> n.right >> n.feature >> n.threshold >> n.value;
            tree.nodes.push_back(n);
        }
        model.trees.push_back(tree);
    }
    return model;
}
