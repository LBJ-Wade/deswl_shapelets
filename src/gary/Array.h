// 	$Id: Array.h,v 1.4 2002/09/10 20:57:51 garyb Exp $
// 2-dimensional array of arbitrary type.  No mathematical intent.
//---------------------------------------------------------------------------
#ifndef ArrayH
#define ArrayH
//---------------------------------------------------------------------------

#include <vector>
#include <assert.h>

namespace laguerre {

    template <class T>
    class Array {

    public:

        explicit Array(uint _n1=0,uint _n2=0) : n1(_n1),n2(_n2),
        data(_n1,std::vector<T>(_n2)) {}
        Array(uint _n1,uint _n2,const T& value) : n1(_n1),n2(_n2),
        data(_n1,std::vector<T>(_n2,value)) {}
        ~Array() {}
        Array(const Array& rhs) : n1(rhs.n1),n2(rhs.n2),data(rhs.data) {}
        Array<T>& operator=(const Array<T>& rhs)
        { 
            if (this == &rhs) return *this; n1 = rhs.n1; n2 = rhs.n2; 
            data = rhs.data; return *this; 
        }
        std::vector<T>& operator[](uint i) {return data[i]; }
        const std::vector<T>& operator[](uint i) const {return data[i];}
        //operator T() const {assert(getN1()==1&&getN2()==1); return data[0][0];}
        uint getN1() const {assert(data.size()==n1); return n1;}
        uint getN2() const 
        { 
            assert(data.size()==0 || data[0].size()==n2);
            return n2;
        }
        const T& get(uint i,uint j) const 
        {assert(i<getN1()); assert(j<getN2()); return data[i][j]; }
        T& get(uint i,uint j) 
        {assert(i<getN1()); assert(j<getN2()); return data[i][j]; }
        const T& operator()(uint i,uint j) const {return get(i,j);}
        T& operator()(uint i,uint j) {return get(i,j);}

        const std::vector<T>& getRow(uint i) const
        {assert(i<getN1()); return data[i]; }
        std::vector<T> getCol(uint j) const
        {
            assert(j<getN2()); 
            std::vector<T> temp(getN1()); 
            for(uint i=0;i<getN1();i++) temp[i] = data[i][j]; 
            return temp; 
        }
        void SetN1(uint newn1);
        void SetN2(uint newn2);
        void SetAllValues(const T& value)
        {
            for(uint i=0;i<getN1();i++) for(uint j=0;j<getN2();j++) 
                data[i][j] = value; 
        }
        void Set(uint i,uint j,const T& value) {data[i][j] = value;}
        void SetRow(uint i,const std::vector<T>& row) 
        {
            assert(i<getN1());
            assert(row.size()==getN2());
            for(uint j=0;j<getN2();j++) data[i][j] = row[j]; 
        }
        void SetCol(uint j,const std::vector<T>& col)
        {
            assert(j<getN2());
            assert(col.size()==getN1());
            for(uint i=0;i<getN1();i++) data[i][j] = col[i]; 
        }
        void SetRow(uint i,const std::vector<T>& row,uint start,uint end) 
        { 
            assert(i<getN1());
            assert(start <= end);
            assert(end < row.size());
            assert(end-start+1==getN2());
            for(uint j=0;j<getN2();j++) data[i][j] = row[j+start]; 
        }
        void SetCol(uint j,const std::vector<T>& col,uint start,uint end)
        {
            assert(j<getN2());
            assert(start <= end);
            assert(end < col.size());
            assert(end-start+1==getN1());
            for(uint i=0;i<getN1();i++) data[i][j] = col[i+start]; 
        }
        void AddRow(std::vector<T> row) 
        {
            assert(row.size()==getN2());
            data.push_back(row); n1++; 
        }
        void AddCol(std::vector<T> col) 
        {
            assert(col.size()==getN1());
            for(uint i=0;i<getN1();i++) data[i].push_back(col[i]); n2++; 
        }
        void Reserve(uint n1r,uint n2r);

    protected :

        uint n1,n2;
        std::vector<std::vector<T> > data;
    };

    template <class T>
    inline void Array<T>::SetN1(uint newn1)
    {
        if (n1 == 0) 
            data = std::vector<std::vector<T> >(newn1,std::vector<T>(n2));
        else if (newn1 > n1) 
            data.insert(data.end(),newn1-n1,std::vector<T>(n2));
        else data.erase(data.begin()+newn1,data.end());
        n1 = newn1;
    }

    template <class T>
    inline void Array<T>::SetN2(uint newn2)
    {
        if (n2 == 0) for(uint i=0;i<n1;i++) data[i] = std::vector<T>(newn2);
        else if (newn2 > n2) for(uint i=0;i<n1;i++)
            data[i].insert(data[i].end(),newn2-n2,0.);
        else for(uint i=0;i<n1;i++)
            data[i].erase(data[i].begin()+newn2,data[i].end());
        n2 = newn2;
    }

    template <class T>
    inline void Array<T>::Reserve(uint n1r,uint n2r)
    {
        if (n1r > n1) data.reserve(n1r);
        if (n2r > n2) for(uint i=0;i<n1;i++) data[i].reserve(n2r);
    }

}

#endif
