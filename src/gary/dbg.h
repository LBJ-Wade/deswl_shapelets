//---------------------------------------------------------------------------
// Debugging and exception classes:

namespace laguerre {

    inline void error(const char *s,const char *s2 = "")
    {
        cerr << "Error: " << s << ' ' << s2 << endl;
        exit(1);
    }

    // Define a base exception class with error message:
    class MyException 
    {
    public:
        MyException(const string &m=""): msg(m) {};
        MyException(const char* c): msg(c) {};
        void set(const string &m="") {msg=m;}
        void set(const char* c) {msg=c;}
        void dump(ostream &os) const {os << "Error: " << msg << endl;}
        void quit(const int exit_code=1) {
            cerr << "Error: " << msg << endl;
            exit(exit_code);
        }
        string msg;
    };

}
#endif
