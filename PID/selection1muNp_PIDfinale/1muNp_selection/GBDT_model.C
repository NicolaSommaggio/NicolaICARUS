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

int main() {
    GBDTModel model = load_model("GBDT_small_stat_newTS_newvars_onlycoll.txt");
    cout << "learning_rate = " << model.learning_rate << endl;
    cout << "n_classes = " << model.n_classes << endl;
    cout << "n_trees = " << model.trees.size() << endl;

    //ifstream test("y_pred.txt");
    //std::string line;
    //ofstream check_pred("check_pred.txt");
    //ofstream check_pred_py("check_pred_py.txt");
    //while(std::getline(test,line))
      //{
      //double lr1, lr2, lr3, lr4, lr5, lr6, lr7, lr8, lr9, lr10, lr11, lr12, lr13, lr14, lr15, depE, depE_daughter, angle_daughter, pred;
      //std::stringstream ss(line);
      //ss >> lr1 >> lr2 >> lr3 >> lr4 >> lr5 >> lr6 >> lr7 >> lr8 >> lr9 >> lr10 >> lr11 >> lr12 >> lr13 >> lr14 >> lr15 >> depE >> depE_daughter >> angle_daughter >> pred;
      //std::vector<double> x = {lr1, lr2, lr3, lr4, lr5, lr6, lr7, lr8, lr9, lr10, lr11, lr12, lr13, lr14, lr15, depE, depE_daughter, angle_daughter, pred};
      //check_pred << model.predict_class(x) << endl;
      //check_pred_py << pred << endl; 
      //}

    //std::vector<double> x = {0.247, 0.09, 0.2, 0.009, 0.163, -0.166, -0.053, -0.24, -0.094, 0.115, -0.081, 0.076, -0.192, -0.041, 0.155, 3.3612685799598694, -1, -1};
    std::vector<double> x = {0.247, 0.09, 0.2, 0.009, 0.163, -0.166, 0.053, -0.24, -0.094, 0.115, -0.081, 0.076, -0.192, -0.041, 0.155, 0.3612685799598694, -1, -1};
    
    
    auto proba = model.predict_proba(x);
    for (int k = 0; k < model.n_classes; k++)
    {
    cout << fixed << setprecision(8);
    cout << "P(classe " << k << ") = " << proba[k] << endl;
    }
    cout << "Classe predetta: " << model.predict_class(x) << endl;
}
