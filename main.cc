#include "structures.h"
#include "obj_parser.h"



int main(int argc, char *argv[]) {

    int fCount = argc;

    std::vector<Mesh> Models;

    while (--fCount) {
        Models.push_back(getModel(argv[fCount]));
    }

}
