#ifndef __IDIGM_FILEIO_H__
#define __IDIGM_FILEIO_H__

#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include "util.cpp"
#include "montgomery.cpp"
#include "factorgraph.cpp"
#include "crt.cpp"


str_t removeLeading(vec_t<char> charsToRemove, const str_t& s)
{
    uint_t n = s.size();
    uint_t resultStartPos {0};
    for (resultStartPos = 0; resultStartPos < n; ++resultStartPos) {
        if ( ! isIn(charsToRemove, s[resultStartPos]))
            break;
    }
    return str_t(s, resultStartPos);
}

str_t removeLeadingWhitespace(const str_t& s)
{
    return removeLeading({' ', '\n', '\t'}, s);
}

str_t takeUntil(vec_t<char> separatorSyms, const str_t& str)
{
    uint_t n = str.size();
    uint_t endPos {0};
    for (endPos = 0; endPos < n; ++endPos) {
        if (isIn(separatorSyms, str[endPos]))
            break;
    }
    return str.substr(0, endPos);
}

str_t getNthWord(str_t str, uint_t n=0,
                 vec_t<char> separatorSyms={' ', '\n', '\t', ',', ';'})
{
    str_t word("");
    for (uint_t i = 0; i <= n; ++i) {
        str = removeLeading(separatorSyms, str);
        word = takeUntil(separatorSyms, str);
        if (word == "")
            break;
        str = str_t(str, word.size());
    }
    return word;
}

uint_t getWordCount(str_t s)
{
    uint_t wordCount {0};
    while (getNthWord(s, wordCount) != "")
        wordCount += 1;
    return wordCount;
}

bool fileExists(str_t filename)
{
    struct stat buf;
    return (stat(filename.c_str(), &buf) == 0);
}


class FileHandler {
private:
    FILE* file;
    char* buffer;
    str_t bufferStr;
    str_t currentToken;
    static const uint_t bufferSize = 1U<<27;
    static const uint_t bufferStrReloadThreshold = 1U<<10;
    
    void reloadBufferStr()
    {
        if (feof(file))
            return;
        int numBytesRead = fread(buffer, sizeof(char), bufferSize, file);
        if (numBytesRead < 1)
            assert (!ferror(file)); //TODO: handle error instead of aborting
        buffer[numBytesRead] = '\0';
        bufferStr += str_t(buffer);
    }

public:
    FileHandler(str_t filename) :
        file{nullptr}, buffer{nullptr}, bufferStr{""}, currentToken{""}
    {
        if ( ! fileExists(filename))
            file = fopen(filename.c_str(), "w+");
        else
            file = fopen(filename.c_str(), "r+");
        assert ( file != nullptr );
        buffer = new char[bufferSize+1];
        buffer[bufferSize] = '\0';
        nextToken();
    }
    
    FileHandler(FILE* f) :
        file{f}, buffer{nullptr}, bufferStr{""}, currentToken{""}
    {
        assert ( file != nullptr );
        assert ( fseek(file, 0, SEEK_SET) == 0 );
        buffer = new char[bufferSize+1];
        buffer[bufferSize] = '\0';
        nextToken();
    }
    
    FileHandler(const FileHandler& other) = delete;
    FileHandler& operator=(const FileHandler& other) = delete;
    
    ~FileHandler()
    {
        if (file != nullptr) fclose(file);
        delete[] buffer;
    }
    
    //NOTE: Caller's responsibility to ensure that no file operations are
    // performed after calling close.
    void close() {
        fclose(file);
        file = nullptr;
        bufferStr = "";
        currentToken = "";
    }
    
    bool atEnd() {
        return feof(file)
               && (bufferStr.size() == 0)
               && (currentToken.size() == 0);
    }
    
    str_t getCurrentToken() { return currentToken; }
    
    // position (in bytes from file beginning) of first char of currentToken
    long getCurrentPosition() {
        return ftell(file) - (long)(bufferStr + currentToken).size();
    }
    
    void reset() {
        assert ( fseek(file, 0, SEEK_SET) == 0 );
        bufferStr = "";
    }
    
    void toBeginning() {
        reset();
        nextToken();
    }
    
    void toEnd() {
        assert ( fseek(file, 0, SEEK_END) == 0 );
        currentToken = "";
        bufferStr = "";
    }
    
    str_t nextToken() {
        if (bufferStr.size() < bufferStrReloadThreshold)
            reloadBufferStr();
        bufferStr = removeLeadingWhitespace(bufferStr);
        char wordBuffer[1U<<12];
        int numArgsRead = sscanf(bufferStr.c_str(), "%s", wordBuffer);
        if (numArgsRead < 1) {
            currentToken = "";
            return "";
        }
        currentToken = str_t(wordBuffer);
        bufferStr = str_t(bufferStr, currentToken.size());
        return currentToken;
    }
    
    bool toNextOccurrence(str_t target) {
        long triesLeft {1L<<34};
        while ( !atEnd() && (currentToken != target)) {
            nextToken();
            assert ( --triesLeft >= 0 );
        }
        return currentToken == target;
    }
    
    long getNextOccurrencePos(const str_t& target) {
        if ( ! toNextOccurrence(target))
            return -1;
        return getCurrentPosition();
    }
    
    long getFirstOccurrencePos(const str_t& target) {
        long originalPos = ftell(file);
        assert ( originalPos >= 0 );
        str_t originalBufferStr = bufferStr;
        str_t originalCurrentToken = currentToken;
        
        toBeginning();
        long targetFirstPos = getNextOccurrencePos(target);
        
        assert ( fseek(file, originalPos, SEEK_SET) == 0 );
        bufferStr = originalBufferStr;
        currentToken = originalCurrentToken;
        
        return targetFirstPos;
    }
    
    // (over-)write at beginning of currentToken
    void writeAtCurrentPos(str_t str) {
        long pos = getCurrentPosition();
        assert ( fseek(file, pos, SEEK_SET) == 0 );
        assert ( fwrite(str.c_str(), str.size(), 1, file) == 1 );
        assert ( fflush(file) == 0 );
        bufferStr = "";
        reloadBufferStr();
        currentToken = "";
    }
    
    // overwrite range with whitespace
    void eraseRange(long start, long end) {
        assert ( fseek(file, start, SEEK_SET) == 0 );
        assert ( end >= start );
        long_uint_t n = (long_uint_t)(end - start);
        str_t blankStr(n, ' ');
        assert ( fwrite(blankStr.c_str(), n, 1, file) == 1 );
        assert ( fflush(file) == 0 );
        bufferStr = "";
        reloadBufferStr();
        currentToken = "";
    }
};

bool representsZero(const str_t& s, bool allowSignPrefix=true)
{
    if (getWordCount(s) != 1)
        return false;
    str_t r = getNthWord(s, 0);
    if (r.size() < 1)
        return false;
    if (r[0] == '-' || r[0] == '+') {
        if ( ! allowSignPrefix)
            return false;
        return representsZero(str_t(r, 1), false);
    }
    if (r[0] != '0')
        return false;
    r = removeLeading({'0'}, r);  // NOTE: allow multiple leading "0"s
    if (r.size() == 0)
        return true;
    if (r.size() == 1 || r[0] != '.')  // NOTE: disallow ".0" and "0."
        return false;
    r = str_t(r, 1);
    r = removeLeading({'0'}, r);
    return (r.size() == 0);
}

long atolErr(const str_t& s)
{
    assert ( s.size() >= 1 );
    long result = std::atol(s.c_str());
    // atol returns 0 both when s represents zero and when s is not convertible
    if (result == 0)
        assert ( representsZero(s) );
    return result;
}

float atofErr(const str_t& s)
{
    assert ( s.size() >= 1 );
    float result = std::atof(s.c_str());
    // atof returns 0 both when s represents zero and when s is not convertible
    if (result == 0.0)
        assert ( representsZero(s) );
    return result;
}

// Declaration without definition, to be specialized for types of interest
template <typename T>
T deserialize(const str_t& s);

template <>
long deserialize(const str_t& s) { return atolErr(s); }

template <>
int deserialize(const str_t& s)
{
    long x = atolErr(s);
    assert ( x < 1L<<31 );
    assert ( x > -(1L<<31) );
    return (int)x;
}

template <>
float deserialize(const str_t& s) { return atofErr(s); }

template <>
bool deserialize(const str_t& s) { return (s == "true" || s == "1"); }

template <>
long_uint_t deserialize(const str_t& s)
{
    long v = atolErr(s);
    assert ( v >= 0 );
    return (long_uint_t)v;
}

template <>
uint_t deserialize(const str_t& s)
{
    long x = deserialize<long>(s);
    assert ( x >= 0 );
    assert ( x < 1L<<32 );
    return (uint_t)x;
}

template<>
MPInt deserialize(const str_t& s)
{
    mpz_t x; mpz_init_set_str(x, s.c_str(), 10);
    MPInt result(x);
    mpz_clear(x);
    return result;
}

template<>
str_t deserialize(const str_t& s) { return s; }


bool loadObject(FileHandler& fh, str_t& result)
{
    result = fh.getCurrentToken();
    return true;
}

bool loadObject(FileHandler& fh, bool& result)
{
    result = deserialize<bool>(fh.getCurrentToken());
    return true;
}

bool loadObject(FileHandler& fh, long& result)
{
    result = deserialize<long>(fh.getCurrentToken());
    return true;
}

bool loadObject(FileHandler& fh, int& result)
{
    result = deserialize<int>(fh.getCurrentToken());
    return true;
}

bool loadObject(FileHandler& fh, long_uint_t& result)
{
    result = deserialize<long_uint_t>(fh.getCurrentToken());
    return true;
}

bool loadObject(FileHandler& fh, uint_t& result)
{
    result = deserialize<uint_t>(fh.getCurrentToken());
    return true;
}

bool loadObject(FileHandler& fh, float& result)
{
    result = deserialize<float>(fh.getCurrentToken());
    return true;
}

bool loadObject(FileHandler& fh, Var& result)
{
    if (fh.getCurrentToken() != "Var")
        // TODO: Log error: could not find a Var at current location
        return false;
    label_t label = deserialize<long_uint_t>(fh.nextToken());
    long_uint_t size = deserialize<long_uint_t>(fh.nextToken());
    assert ( size < 1UL<<32 );
    result = Var(label, (uint_t)size);
    return true;
}

bool loadObject(FileHandler& fh, MPInt& result)
{
    str_t s = fh.getCurrentToken();
    result = deserialize<MPInt>(fh.getCurrentToken());
    return true;
}

template <typename P>
bool loadObject(FileHandler& fh, Zp<P>& result)
{
    long rawValue = deserialize<long>(fh.getCurrentToken());
    assert ( rawValue >= 0 ); assert ( rawValue < 1L<<32 );
    // NOTE: Assuming Zp values are stored "raw", i.e. in Montgomery form
    result = Zp<P>((scalar_t)rawValue, true);
    return true;
}

template <typename T>
bool loadObject(FileHandler& fh, vec_t<T>& result)
{
    if (fh.getCurrentToken() != "{") {
        if (fh.nextToken() != "{")
            // TODO: Log error: could not find a vector at current location
            return false;
    }
    T x;
    while ( ! fh.atEnd()) {
        str_t token = fh.nextToken();
        if (token == "}") 
            break;
        if ( ! loadObject(fh, x))
            // TODO: Log error: could not load element
            return false;
        result.push_back(x);
    }
    return true;
}

template <typename T>
bool loadObject(FileHandler& fh, set_t<T>& result)
{
    vec_t<T> vecResult;
    bool success = loadObject(fh, vecResult);
    result = vectorToSet(vecResult);
    return success;
}

template <typename K, typename V>
bool loadObject(FileHandler& fh, std::pair<K,V>& result)
{
    if (fh.getCurrentToken() != "{") {
        if (fh.nextToken() != "{")
            // TODO: Log error: could not find a pair at current location
            return false;
    }
    
    fh.nextToken();
    K first;
    if ( ! loadObject(fh, first))
        return false;
    
    if (fh.nextToken() == ",")
        fh.nextToken();
    
    V second;
    if ( ! loadObject(fh, second))
        return false;
    
    if ( ! (fh.nextToken() == "}"))
        return false;
    
    result = std::pair(first, second);
    return true;
}

template <typename K, typename V>
bool loadObject(FileHandler& fh, map_t<K,V>& result)
{
    vec_t<std::pair<K,V>> vecOfPairs;
    if ( ! loadObject(fh, vecOfPairs))
        return false;
    for (std::pair<K,V> p : vecOfPairs)
        result.insert(p);
    return true;
}

template <typename T>
bool loadObject(FileHandler& fh, Factor<T>& result)
{
    if ( ! fh.toNextOccurrence("Factor") )
        // TODO: Log error: could not find a Factor
        return false;
    
    if (fh.nextToken() != "label")
        // TODO: Log error: could not find Factor label field
        return false;
    label_t label = deserialize<label_t>(fh.nextToken());
    
    if (fh.nextToken() != "vars")
        // TODO: Log error: could not find Factor vars field
        return false;
    vec_t<Var> vars;
    if ( ! loadObject(fh, vars))
        // TODO: Log error: could not load vars
        return false;
    
    if (fh.nextToken() != "data")
        // TODO: Log error: could not find Factor data
        return false;
    vec_t<T> data;
    if ( ! loadObject(fh, data))
        // TODO: Log error: could not load Factor data
        return false;
    
    result = Factor<T>(label, vars, data.data());
    return true;
}

template <typename T>
bool loadObject(FileHandler& fh, FactorGraph<T>& result)
{
    if (fh.getCurrentToken() != "FactorGraph") {
        if (fh.nextToken() != "FactorGraph")
            // TODO: Log error: could not find a FactorGraph
            return false;
    }
    
    if (fh.nextToken() != "boundary")
        // TODO: Log error: could not find boundary field
        return false;
    set_t<Var> boundary;
    if ( ! loadObject(fh, boundary))
        // TODO: Log error: could not load boundary
        return false;
    
    if (fh.nextToken() != "factors")
        // TODO: Log error: could not find list/vector of factors
        return false;
    set_t<Factor<T>> factors;
    if ( ! loadObject(fh, factors))
        // TODO: Log error: could not load factors
        return false;
    
    result = FactorGraph<T>(factors, boundary);
    return true;
}

template <typename T>
bool loadFactorGraph(FileHandler& fh, FactorGraph<T>& result)
{
    if ( ! fh.toNextOccurrence("FactorGraph"))
        // TODO: Log error: could not find a FactorGraph
        return false;
    
    return loadObject(fh, result);
}

// Load only the boundary of the next FactorGraph
bool loadBoundary(FileHandler& fh, vec_t<Var>& boundary)
{
    if ( ! fh.toNextOccurrence("FactorGraph"))
        // TODO: Log error: could not find a FactorGraph
        return false;
    if (fh.nextToken() != "boundary")
        // TODO: Log error: could not find boundary field
        return false;
    if ( ! loadObject(fh, boundary))
        // TODO: Log error: could not load boundary
        return false;
    return true;
}


template <typename T>
str_t serialize(const T& x)
{
    return std::to_string(x);
}

str_t serialize(const float& f)
{
    char buf[64];
    sprintf(buf, "%.32f", f);
    return str_t(buf);
}

str_t serialize(const double& f)
{
    char buf[128];
    sprintf(buf, "%.64f", f);
    return str_t(buf);
}

str_t serialize(const Var& v)
{
    return "Var " + serialize(v.label) + " " + serialize(v.size);
}

template <typename P>
str_t serialize(const Zp<P>& x)
{
    // NOTE: Assuming Zp values are stored "raw", i.e. in Montgomery form
    return serialize<scalar_t>(x.raw());
}

str_t serialize(const MPInt& x) { return MPInt::to_string(x); }

str_t serialize(const str_t& s) { return s; }

void writeObject(FileHandler& fh, const int& x)
{ fh.writeAtCurrentPos(serialize(x)); }

void writeObject(FileHandler& fh, const long& x)
{ fh.writeAtCurrentPos(serialize(x)); }

void writeObject(FileHandler& fh, const float& x)
{ fh.writeAtCurrentPos(serialize(x)); }

void writeObject(FileHandler& fh, const uint_t& x)
{ fh.writeAtCurrentPos(serialize(x)); }

void writeObject(FileHandler& fh, const long_uint_t& x)
{ fh.writeAtCurrentPos(serialize(x)); }

void writeObject(FileHandler& fh, const str_t& x)
{ fh.writeAtCurrentPos(x); }

void writeObject(FileHandler& fh, const MPInt& x)
{ fh.writeAtCurrentPos(serialize(x)); }

void writeObject(FileHandler& fh, const Var& x)
{ fh.writeAtCurrentPos(serialize(x)); }

template <typename P>
void writeObject(FileHandler& fh, const Zp<P>& x)
{ fh.writeAtCurrentPos(serialize(x)); }

template <typename K, typename V>
void writeObject(FileHandler& fh, const std::pair<K,V>& p)
{
    fh.writeAtCurrentPos("{ ");
    writeObject(fh, p.first);
    fh.writeAtCurrentPos(" ");
    writeObject(fh, p.second);
    fh.writeAtCurrentPos(" } ");
}

template <typename T>
void writeObject(FileHandler& fh, const vec_t<T>& v)
{
    fh.writeAtCurrentPos("{ ");
    for (T x : v) {
        writeObject(fh, x);
        fh.writeAtCurrentPos(" ");
    }
    fh.writeAtCurrentPos("} ");
}

template <typename K, typename V>
void writeObject(FileHandler& fh, const map_t<K,V>& m)
{
    fh.writeAtCurrentPos("{ ");
    for (std::pair<K,V> p : m)
        writeObject(fh, p);
    fh.writeAtCurrentPos(" } ");
}

/* If atCurrentLocation, (over-)write starting at current fh position;
 * else write to end of file. */
template <typename F>
void writeFactorToFile(FileHandler& fh, const Factor<F>& factor,
                       bool atCurrentLocation=false, uint_t indent=0)
{
    if ( ! atCurrentLocation) {
        fh.toEnd();
        fh.writeAtCurrentPos("\n");
    }
    str_t indentStr;
    for (uint_t i = 0; i < indent; ++i)
        indentStr += " ";
    fh.writeAtCurrentPos(indentStr + "Factor\n");
    fh.writeAtCurrentPos(indentStr +"  label "+ serialize(factor.label) +"\n");
    fh.writeAtCurrentPos(indentStr +"  vars {\n"+ indentStr);
    for (Var v : factor.getVars())
        fh.writeAtCurrentPos("    " + serialize(v));
    fh.writeAtCurrentPos("\n" + indentStr + "  }\n");
    fh.writeAtCurrentPos(indentStr +"  data {\n    "+ indentStr);
    for (long_uint_t k = 0; k < factor.getVolume(); ++k)
        fh.writeAtCurrentPos(serialize(factor.valueAt(k)) + " ");
    fh.writeAtCurrentPos("\n" + indentStr + "  }\n");
}

template <typename F>
void writeFGToFile(FileHandler& fh, const FactorGraph<F>& g,
                   bool atCurrentLocation=false)
{
    if ( ! atCurrentLocation) {
        fh.toEnd();
        fh.writeAtCurrentPos("\n");
    }
    fh.writeAtCurrentPos("FactorGraph\n");
    fh.writeAtCurrentPos("  boundary {\n");
    for (Var v : g.getBoundary())
        fh.writeAtCurrentPos("    " + serialize(v));
    fh.writeAtCurrentPos("\n  }\n");
    fh.writeAtCurrentPos("  factors {\n");
    for (Factor<F> f : g.getFactors())
        writeFactorToFile(fh, f, true, 4);
    fh.writeAtCurrentPos("\n  }\n");
}

template <typename T>
void writeParamToFile(FileHandler& fh,
                      const str_t& paramName, const T& paramValue)
{
    long preexistingPos = fh.getFirstOccurrencePos(paramName);
    // if param is already defined then erase the old definition
    if (preexistingPos >= 0) {
        // TODO: Log warning: overwriting exisiting parameter value
        // go to beginning of old param value
        fh.toBeginning(); fh.toNextOccurrence(paramName); fh.nextToken();
        T oldValue;
        // read one token past old value
        assert ( loadObject(fh, oldValue) ); fh.nextToken();
        long oldValueEnd = fh.getCurrentPosition() - 1;
        fh.eraseRange(preexistingPos, oldValueEnd);
    }
    
    // append new param definition to end of file
    fh.toEnd();
    fh.writeAtCurrentPos("\n" + paramName + " ");
    writeObject(fh, paramValue);
    fh.writeAtCurrentPos("\n");
}

template <typename T>
bool loadParamFromFile(FileHandler& fh,
                       const str_t& paramName, T& result)
{
    fh.toBeginning();
    if ( ! fh.toNextOccurrence(paramName))
        return false;
    fh.nextToken();
    return loadObject(fh, result);
}

template <typename T>
bool copyParam(FileHandler& src, FileHandler& dst, str_t paramName)
{
    T x;
    if (loadParamFromFile(src, paramName, x)) {
        writeParamToFile(dst, paramName, x);
        return true;
    }
    return false;
}


#endif  // __IDIGM_FILEIO_H__

