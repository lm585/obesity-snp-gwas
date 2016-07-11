/*
input file 1 (vcf file)
$ head  obes.rsID.genotypes.vcf   | cut -f 1-7,9-20
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  FORMAT  HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00105 HG00106 HG00107 HG00108
10      4655565 rs10458787      A       G       100     PASS    GT      1|0     1|0     1|1     1|1     1|1     1|0     0|1     0|0     1|0     1|0     1|1
10      4920413 rs10904363      G       C       100     PASS    GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0
10      16299951        rs10508503      C       T       100     PASS    GT      0|0     0|0     0|0     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0
10      23858211        rs16923476      G       A       100     PASS    GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0
10      31990623        rs7081678       G       A       100     PASS    GT      0|0     1|0     0|0     1|0     0|0     0|0     0|0     1|0     0|0     0|0     0|0
10      45133277        rs11239187      A       G       100     PASS    GT      1|0     1|1     0|0     0|0     0|1     0|0     0|1     1|0     0|1     1|1     0|1
10      78646536        rs2116830       G       T       100     PASS    GT      0|1     0|0     1|0     1|0     0|1     0|0     0|1     1|1     0|0     0|0     0|0
10      104906211       rs11191580      T       C       100     PASS    GT      0|0     1|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0     0|0     1|0
11      2839751 rs2237892       C       T       100     PASS    GT      0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0

input file 2:
HG00096 
HG00097 
HG00099 
HG00108

input file 3: snp-ID obes-increase-allele (GWAS p-value < 5E-8)
rs10458787	A
rs10904363	C
rs11191580	C
rs2237892	T

output (cumulative obes-increase-alleles):
HG00096	8 8 (4 SNPS, all homozyg risk alleles)
HG00097 2 8 (1 snp with homozyg risk allele, or 2 snps with each having single risk allele)
HG00099 5 8 ()
HG00108 5  6 (./. no call for one of the SNPs) 

map<string, string> refIDmap; //"rs10458787" -> "AG" ; [0] ref base; [1] 1st alt base
vector<people> peopleVect; //
vector<string> inputSNPsVect;
vector<char> inputSNPsRiskAllele;
map<string, int> inputPeopleMap; //"HG00096" -> 0; "HG00097" -> 0

for each person in the peopleVect contained in the inputPeopleMap
  the person is in VCF file
  c1 = 0; c2 = 0;
  for each snp in the inputSNPsVect
    if snp NOT found in the person's record, exit;
    if snp found
      if !(inputSNPsRiskAllele being one of two alleles in 1000-genome)
         output snp inputSNPsRiskAllele not in the 1000-genome
         exit
      1|0, 1|1, 0|0, 0|1, i/j, ./. 
      if . in genotype
        doing nothing
      else if inputSNPsRiskAllele is ref allele (0)
        c2 += 2; 
        count how many '0' (0/1/2)
        c1 += how-many
      else if inputSNPsRiskAllele is ref allele (1)
        c2 += 2;
        count how many '1' (0/1/2)
        c1 += how-many
      else;
  output personID	c1	c2

output if any person NOT in VCF file
 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <sstream>
#include <map>

using namespace std;

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);


struct people {
 string id; //HG00096
 map<string, string> genoMap; //"rs10458787" -> 1|0; "rs10904363" -> 0|0
};

int main(int argc, char *argv[])
{
 if(argc != 6)
 {
  cerr << argv[0] << "   1[VCF file]  2[a list of individual ID (HG00096)]     3[obes rs ID file with obes increase allele with GWAS p < 5E-8 (rs1011731	A)]       4[# of fields in header/data line]  5[output-file]" << endl;
  cerr << "******************************************************************************" << endl;
  cerr<<"warning: all VCF files are assumed to have same format" << endl;
  cerr << "Format (order and number of fields are assumed): " << endl;
  cerr<<"#CHROM POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096 HG00097 HG00099 HG00100 HG00101 ..." << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;

 map<string, int> inputPeopleMap; //"HG00096" -> 0; "HG00097" -> 0
 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent2(lineFields, " \t\n", line);
  for(int i = 0; i < lineFields.size(); i++)
  {
   inputPeopleMap[lineFields[i] ] = 0;
  }
  getline(input, line);
 }
 input.close();

 vector<string> inputSNPsVect;
 vector<char> inputSNPsRiskAllele;
 readFile(argv[3], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent2(lineFields, " \t\n", line);
  if(lineFields.size() != 2)
  {
   cerr << line << " input SNP line format wrong. \'rs10904363	A\' " << endl;
   exit(1);
  }
  inputSNPsVect.push_back(lineFields[0]);
  inputSNPsRiskAllele.push_back(lineFields[1][0]);
  getline(input, line);
 }
 input.close();

  
 writeFile(argv[5], output);
 map<string, string> refIDmap; //"rs10458787" -> "AG" ; [0] ref base; [1] 1st alt base
 vector<people> peopleVect; //
 vector<string> vcfHeader;
 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  if(lineFields.size() == atoi(argv[4])  && lineFields[0] == "#CHROM" )
  {
   for(int i = 0; i < lineFields.size(); i++)
   {
    vcfHeader.push_back(lineFields[i]);
   }
   for(int i = 9; i < lineFields.size(); i++)
   {
    people peop;
    peop.id = lineFields[i];
    peopleVect.push_back(peop);
   }
  }
//0       1        2             3        4      5         6      7    8       9
//#CHROM  POS     ID             REF     ALT     QUAL    FILTER  info  FORMAT  HG00096
//10      4655565 rs10458787      A       G       100     PASS   xxxx  GT      1|0
  else if (lineFields.size() == atoi(argv[4])) //body line
  {
   string geno;
   people peop;
   geno = "";
   geno += lineFields[3][0];
   geno += lineFields[4][0];
   refIDmap[lineFields[2] ] = geno;
   for(int i = 9; i < lineFields.size(); i++)
   {
    peopleVect[i-9].genoMap[lineFields[2] ] = lineFields[i]; //peopleVect[0] id: lineFields[9];  "rs10458787" -> 1|0
                                                             //          [1] id: lineFields[10]; "rs10458787" -> x|y   
   }
  }
  else;
  getline(input, line);
 }
 input.close();

 for(int i = 0; i < peopleVect.size(); i++)
 {
  if(inputPeopleMap.count(peopleVect[i].id) != 0)
  {
   inputPeopleMap[peopleVect[i].id ] = 1;
   string  genoStr;
   int c1, c2;
   string snp;
   char inputAllele;
   vector<string> al12vect;

   c1 = 0;
   c2 = 0;
   for(int j = 0; j < inputSNPsVect.size(); j++)
   {
    snp = inputSNPsVect[j];
    inputAllele = inputSNPsRiskAllele[j];
    if(peopleVect[i].genoMap.count(snp) != 0)
    {
     //snp also contained in refIDmap
     if(refIDmap[snp].find(inputAllele) == std::string::npos)
     {
      cerr << snp << " inputSNPsRiskAllele " << inputAllele << " not in the 1000-genome" << endl;
      exit(1);
     }
     genoStr = peopleVect[i].genoMap[snp]; //0|1
     getFieldContent2(al12vect, "|/ \n", genoStr); //0,1
     //input risk allele being one of the two alleles
     if(al12vect.size() >= 2)
     {
      if(al12vect[0] ==  "." || al12vect[1] ==  ".")
      {
       //doing nothing
      }
      else if(inputAllele == refIDmap[snp][0]) //inputSNPsRiskAllele is ref allele (0)
      {
       c2 += 2;
       if(al12vect[0]  == "0")
         c1++;
       if(al12vect[1]  == "0")
         c1++;
      }
      else if(inputAllele == refIDmap[snp][1]) //inputSNPsRiskAllele is ALT allele (1)
      {
       c2 += 2;
       if(al12vect[0]  == "1")
         c1++;
       if(al12vect[1]  == "1")
         c1++;
      }
      else;
     }
    }
    else
    {
     cerr << snp << " not in the VCF file. Exit..." << endl;
     exit(1);
    }
   }
   // if(str1.length() == inputSNPsVect.size() && str2.length() == inputSNPsVect.size())
   {
    output << peopleVect[i].id << '\t' << c1 << '\t' << c2 << endl;
   }
  } //end of if block ------- if(inputPeopleMap.count(peopleVect[i].id) != 0)
 } //end of for loop --- for(int i = 0; i < peopleVect.size(); i++) 
 output.close();
 for(map<string, int>::iterator it = inputPeopleMap.begin(); it != inputPeopleMap.end(); it++)
 {
  if(it->second == 0)
  {
   cerr << it->first << " not in the VCF file" << endl;
  }
 }
 return 0;
}

void readFile(const char * argv, ifstream & input)
{
 input.open(argv, ios::in);
 if(!input)
 {
  cerr << argv << " can not be opened for reading.\n";
  exit(1);
 }
}

void writeFile(const char * argv, ofstream & output)
{
 output.open(argv, ios::out);
 if(!output)
 {
  cerr << argv << " can not be opened for writing.\n";
  exit(1);
 }
}

void getFieldContent(vector<string> & fields, char del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del == line[i])
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   str = line.substr(pos[i] + 1, len);
   fields.push_back(str);
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}

void getFieldContent2(vector<string> & fields, string del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del.find(line[i]) != string::npos)
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   if(len > 0)
   {
    str = line.substr(pos[i] + 1, len);
    fields.push_back(str);
   }
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}


