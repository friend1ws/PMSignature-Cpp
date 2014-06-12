#ifndef __MUTINFO_H__
#define __MUTINFO_H__

struct MutInfo {
	std::string mutationID;
	int sampleID;
	std::string sequence;
	char altBase;
	public:
	MutInfo(std::string mID, int sID, std::string seq, char alt) : mutationID(mID), sampleID(sID), sequence(seq), altBase(alt) {};
	virtual ~MutInfo() {}
};

struct MutInfo_simple {
    int sampleID;
    int patternVec;
    public:
    MutInfo_simple(int sID, int pVec) : sampleID(sID), patternVec(pVec) {};
    virtual ~MutInfo_simple() {}
};


#endif
