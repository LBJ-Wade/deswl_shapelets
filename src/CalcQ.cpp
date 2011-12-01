#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <map>
#include <vector>
#include <cassert>
#include <complex>

// Official code from Sarah:
// tot=0;
// n=0;
// for i=1:nset
//  dg1 = my_g1(i) - true_g1(i);
//  dg2 = my_g2(i) - true_g2(i);
//  tot = tot + dg1^2 + dg2^2;
//  n = n + 2;
// end
// meandgsq=tot/n;
// Q = 1e-4 / meandgsq;


int main(int argc, char* argv[])
{
    if (argc < 3 || argc > 5) {
        std::cerr<<"Usage: calcq correctfile myanswers [nmax] [cheat]\n";
        std::cerr<<"If nmax is specified, only do the first nmax lines\n";
        std::cerr<<"If cheat is specified, use the cheat file to get galaxy types\n";
        exit(1);
    }
    std::ifstream correct(argv[1]);
    std::ifstream myanswers(argv[2]);
    int nmax = (argc >= 4) ? strtol(argv[3],0,0) : -1;
    if (nmax == 0) nmax = -1;
    bool haveCheat = argc >= 5;

    std::ifstream cheat;
    std::vector<int> fileType;
    std::vector<std::string> typeName;
    std::map<std::string,int> typeMap;
    std::vector<int> fileGroup;
    std::vector<int> groupType;
    std::vector<int> groupBatch;
    std::map<std::pair<int,int>,int> groupMap;
    std::vector<std::complex<double> > fileEm;
    std::vector<std::complex<double> > fileEc;
    std::vector<std::string> fileId;
    std::map<std::string,int> idMap;
    int nType = 0;
    int nGroup = 0;
    int nFile = 0;

    if (haveCheat) {
        cheat.open(argv[4]);
        std::string id, doc;
        double e1_c,e2_c;
        while (cheat >> id >> e1_c >> e2_c >> doc) {
            //std::cout<<id<<"  "<<doc<<std::endl;
            size_t k0 = id.find("set");
            k0 += 3;
            id = std::string(id,k0,id.size()-k0);

            size_t k1 = doc.find('_');
            ++k1;
            size_t k2 = doc.find('_',k1);
            std::string name(doc,k1,k2-k1);
            if (typeMap.find(name) == typeMap.end()) {
                typeMap[name] = nType++;
                typeName.push_back(name);
            }
            int type = typeMap[name];
            fileType.push_back(type);
            k1 = doc.find("shear");
            k1 += 6;
            k2 = doc.find('_',k1);
            int batch = strtol(doc.c_str()+k1,0,0);
            //std::cout<<fileType.size()<<" is type "<<type<<" = "<<typeName[type]<<std::endl;
            //std::cout<<fileType.size()<<" is batch "<<batch<<std::endl;
            int k=nFile++;
            idMap[id] = k;
            fileId.push_back(id);

            // For some reason e1 seems to be flpped in these files.
            fileEc.push_back(std::complex<double>(-e1_c,e2_c));
            std::pair<int,int> groupIds(type,batch);
            if (groupMap.find(groupIds) == groupMap.end()) {
                groupMap[groupIds] = nGroup++;
                groupType.push_back(type);
                groupBatch.push_back(batch);
            }
            int group = groupMap[groupIds];
            fileGroup.push_back(group);
        }
        assert(nType == int(typeMap.size()));
        assert(nType == int(typeName.size()));
        assert(nGroup == int(groupMap.size()));
        assert(nGroup == int(groupType.size()));
        assert(nGroup == int(groupBatch.size()));
        assert(nFile == int(idMap.size()));
        assert(nFile == int(fileId.size()));
        assert(nFile == int(fileType.size()));
        assert(nFile == int(fileGroup.size()));
        assert(nFile == int(fileEc.size()));
        std::cout<<"Cheat codes:\n";
        for(int i=0;i<nType;++i) {
            std::cout<<i<<"  =  "<<typeName[i]<<std::endl;
        }
    }

    std::string junk;
    getline(correct,junk);

    std::vector<std::complex<double> > groupDiff(nGroup);
    std::vector<int> nFileInGroup(nGroup);
    double e1_c, e2_c, e1_m, e2_m;
    std::string id1,id2;
    while(myanswers >> id2 >> e1_m >> e2_m >> junk >> junk) {
        correct >> id1 >> e1_c >> e2_c;
        std::cout<<id1<<": "<<e1_c<<" "<<e2_c<<" --- "<<e1_m<<" "<<e2_m;
        assert(id1 == id2);
        fileEm.push_back(std::complex<double>(e1_m,e2_m));
        int i;
        if (haveCheat) {
            i = idMap[id1];
            assert(id1 == fileId[i]);
            assert(norm(std::complex<double>(e1_c,e2_c)-fileEc[i]) < 1.e-5);
        } else {
            i = nFile++;
            fileId.push_back(id1);
            fileEc.push_back(std::complex<double>(e1_c,e2_c));
        }
        if (i == nmax) break;
        std::complex<double> diff = fileEc[i] - fileEm[i];
        //std::cout<<"   "<<diff;
        double normdiff = norm(diff);
        //std::cout<<"   "<<normdiff;
        std::cout<<"   "<<2.e-4/normdiff;
        if (haveCheat) {
            groupDiff[fileGroup[i]] += diff;
            ++nFileInGroup[fileGroup[i]];
            std::cout<<"  "<<fileType[i]<<"  "<<fileGroup[i];
        } else {
            groupDiff.push_back(diff);
            nFileInGroup.push_back(1);
            ++nGroup;
        }
        std::cout<<std::endl;
    }
    if (haveCheat) {
        for(int j=0;j<nGroup;++j) {
            if (nFileInGroup[j] > 0) groupDiff[j] /= nFileInGroup[j];
            std::cout<<"Group "<<j<<" has groupDiff = "<<groupDiff[j];
            std::cout<<"  for "<<nFileInGroup[j]<<" files\n";
        }
    }
    nFile = fileEm.size(); // in case nmax < nFile
    assert(fileEm.size() <= fileEc.size());
    assert(fileEm.size() <= fileId.size());
    assert(fileEm.size() <= fileGroup.size());

    if (haveCheat) {
        std::cout<<"By type:\n";
        for(int i=0;i<nType;++i) {
            //std::cout<<"Type "<<i<<std::endl;
            double sumnormdiff = 0.;
            std::complex<double> sumy=0., sumx=0., sumxsq=0., sumxy=0.;
            int ng=0, nf=0;
            for(int j=0;j<nGroup;++j) if (groupType[j] == i) {
                //std::cout<<"  Using group "<<j<<" with normdiff = "<<norm(groupDiff[j])<<std::endl;
                sumnormdiff += norm(groupDiff[j]);
                ++ng;
            }
            for(int k=0;k<nFile;++k) if (fileType[k] == i) {
                sumy += fileEm[k];
                sumx += fileEc[k];
                sumxsq += fileEc[k]*fileEc[k];
                sumxy += fileEm[k]*fileEc[k];
                ++nf;
            }
            double q = 2.e-4 * ng / sumnormdiff;
            std::cout<<i<<"  "<<typeName[i]<<"    ( groups ";
            //std::cout<<"   q = 2.e-4 * "<<ng<<" / "<<sumnormdiff<<std::endl;
            for(int j=0;j<nGroup;++j) if (groupType[j] == i) {
                std::cout<<j<<" ";
            }
            std::cout<<")\n";
            double nfD = nf;
            std::complex<double> m = 
                (sumxy*nfD - sumx*sumy) / (sumxsq*nfD - sumx*sumx)-1.;
            std::complex<double> c = sumy/nfD - m*sumx/nfD;
            std::cout<<"      q = "<<q<<std::endl;
            std::cout<<"      m = "<<m<<std::endl;
            std::cout<<"      c = "<<c<<std::endl;
        }
    }

    // y = mx + c
    // chisq = Sum (yi - m xi - c)^2
    // d/dm = 2 Sum (yi - m xi - c) (-xi) = 0
    // d/dc = 2 Sum (yi - m xi - c) = 0
    // Sum xi yi = m Sum xi^2 + c Sum xi
    // Sum yi = m Sum xi + c N
    // c = Sum yi / N - m Sum xi / N
    // Sum xi yi / N = m Sum xi^2 / N + Sum xi Sum yi / N^2 - m (Sum xi)^2/N^2
    // m = (N Sum xi yi - Sum xi Sum yi) / (N Sum xi^2 - (Sum xi)^2)
    double sumnormdiff=0.;
    std::complex<double> sumy=0., sumx=0., sumxsq=0., sumxy=0.;
    for(int j=0;j<nGroup;++j) {
        sumnormdiff += norm(groupDiff[j]);
    }
    for(int k=0;k<nFile;++k) {
        sumy += fileEm[k];
        sumx += fileEc[k];
        sumxsq += fileEc[k]*fileEc[k];
        sumxy += fileEm[k]*fileEc[k];
    }
    double q = 2.e-4 * nGroup / sumnormdiff;
    std::cout<<"Overall Q = "<<q<<std::endl;

    double nFileD = nFile;
    std::complex<double> m = 
        (sumxy*nFileD - sumx*sumy) / (sumxsq*nFileD - sumx*sumx)-1.;
    std::complex<double> c = sumy/nFileD - m*sumx/nFileD;
    std::cout<<"m = "<<m<<", c = "<<c<<std::endl;
}
