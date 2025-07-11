#include "helper.h"

namespace help {
    // === Convert double to string and remove trailing zeros ===
    std::string doubleToString(double value) {
        std::ostringstream oss;
        oss << value;
        std::string str = oss.str();
        // Remove trailing zeros and decimal point if necessary
        if (str.find('.') != std::string::npos) {
            str = str.substr(0, str.find_last_not_of('0') + 1);
            if (str.back() == '.') {
                str = str.substr(0, str.size()-1);
            }
        }
        return str;
    }
    
    TTree* getTree(const char* fileName, const char* treeName) {
        // === Lecture du fichier ROOT ===
        TFile *file = TFile::Open(fileName, "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Erreur d'ouverture du fichier" << std::endl;
            return nullptr;
        }

        TTree *tree = nullptr;
        TDirectory *dir = (TDirectory*)file->Get("ana");
        if (!dir) {
            std::cerr << "Répertoire 'ana' introuvable" << std::endl;
            file->Close();
            return nullptr;
        }

        tree = (TTree*)dir->Get(treeName);
        if (!tree) {
            std::cerr << "Arbre 'tree' introuvable dans 'ana'" << std::endl;
            file->Close();
            return nullptr;
        }
        return tree;
    }

    std::vector<double> linspace(double start, double end, int num) {
        std::vector<double> vector;
        vector.reserve(num);
        for (int i=0; i<num; ++i) {
            double value = start+i * (end-start)/(num-1);
            vector.push_back(value);
        }
        return vector;
    }

    void linspace(double start, double end, int num, double* array) {
        for (int i=0; i<num; ++i) {
            array[i] = start+i * (end-start)/(num-1);
        }
    }

    bool ends_with(const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() && str.compare(str.size()-suffix.size(), suffix.size(), suffix) == 0;
    }

    void saveData(const std::string& fileName, double data[nR1][nVolume]) {
        std::cout << "Saving data to " << fileName << std::endl;
        std::ofstream file(fileName, std::ios::out);
        if (!file) {
            std::cerr << "Erreur d'ouverture du fichier" << std::endl;
            return;
        }
        for (int i=0; i<nR1; ++i) {
            for (int j=0; j<nVolume; ++j) {
                file << data[i][j] << " ";
            }
            file << std::endl;
        }
        file.close();
    }

    bool isClose(double x1, double y1, double z1,
        double x2, double y2, double z2,
        double tolerance) {
        return (std::abs(x1 - x2) < tolerance &&
                std::abs(y1 - y2) < tolerance &&
                std::abs(z1 - z2) < tolerance);
        }
}