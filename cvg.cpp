//================================================================================================================
// Name        : cvg.cpp
// Author      : Itai Sharon
// Version     : 1.01, 07/May/2013
// Description : Reports the length, %G+C and coverage for each sequence in the input fasta file. 
//             : Output format is 
//             : <scaf-name>	<size>	<%G+C>	<coverage>	<# of Ns>
//================================================================================================================

#include <string>
#include <set>
#include <map>
#include <iostream>
#include "../common/ReadMappingReader.h"
#include "../common/SeqIORead_fasta.h"
#include "../common/String.h"

using namespace std;

/****************************************************************************************************************/
struct Scaf_stats {
	string	scaf_name;
	size_t	size;
	double	gc;
	size_t	Ns;
	Scaf_stats(string _scaf_name, size_t _size, double _gc, size_t _Ns) : scaf_name(_scaf_name), size(_size), gc(_gc), Ns(_Ns), cvg(0)	
		{}
	Scaf_stats() : scaf_name("N/A"), size(0), gc(0), Ns(0), cvg(0)   
		{}
	const Scaf_stats&	operator += (size_t nbps)	{cvg += nbps; return *this;}
	double			coverage() const		{return double(cvg)/size;}
protected:
	size_t	cvg;
};

/****************************************************************************************************************/
int main(int argc, const char* argv[])
{
	if(argc < 3) {
		cerr << endl << "Usage: " << argv[0] << "[-s <max-num-snps>] <fasta-file> <.sam file1> ... <.sam filen>" << endl << endl;
		return -1;
	}
	string assembly_file;
	size_t max_snps = 2;
	set<string>	sam_files;

	if(!strcmp(argv[1], "-s")) {
		max_snps = atoi(argv[2]);
		assembly_file = argv[3];
		for(int i=4; i<argc; i++)
			sam_files.insert(argv[i]);
	}
	else {
		assembly_file = argv[1];
		for(int i=2; i<argc; i++)
			sam_files.insert(argv[i]);
	}

	map<string, Scaf_stats>	scaf_stats;
	cerr << time_since_epoch() << "Reading assembly file (" << assembly_file << ")" << endl;
	Bio::SeqIORead_fasta<Bio::DNASequence>	in(assembly_file);
	Bio::DNASequence*			seq = NULL;
	while(seq = in.next_seq()) {
		scaf_stats.insert(pair<string, Scaf_stats>(seq->display_id(), Scaf_stats(seq->display_id(), seq->seq().size(), Bio::gc(seq->seq()), Bio::Ns(seq->seq()))));
	}
	cerr << time_since_epoch() << "ok" << endl;

	size_t i;
	for(size_t i=2; (i<argc) && (argv[i][0] == '-'); i++) {
		if(argv[i][1] != 's') {
			max_snps = atoi(argv[++i]);
		}
	}

	for(set<string>::const_iterator it=sam_files.begin(); it!=sam_files.end(); it++) {
		string sam_file = *it;
		cerr << time_since_epoch() << "Reading sam file (" << sam_file << ")" << endl;
		cerr << time_since_epoch();
		print_percent(0, 0);

		Bio::SAMReader	reader(sam_file);
		set<string>	seen_reads;
		size_t		num_read = 0;
		size_t		percent_read = 0;

		// Go over all mappings, collect the first mapping (if multiple exist) for each read
		while(reader.good()) {
			if(((++num_read)%(1000) == 0) && (size_t(reader.percent_entries_read()) > percent_read)) {
				print_percent(percent_read, reader.percent_entries_read());
				percent_read = reader.percent_entries_read();
			}
			
			Bio::ReadMappingPtr mapping(reader.next_mapping());
			if(mapping->read_name() == Bio::SAMReader::end().read_name()) {
				break;
			}
			if(mapping->unmapped() || (mapping->num_snps() > max_snps)) {
				continue;
			}
			// If multiple mappings exist we'll take the first one only
			if(seen_reads.find(mapping->read_name()) != seen_reads.end()) {
				continue;
			}
			seen_reads.insert(mapping->read_name());
			scaf_stats[mapping->ref_name()] += mapping->read_length();
		}
		print_percent(percent_read, 100);
		cerr << endl << time_since_epoch() << "ok" << endl;
	}

	cerr << time_since_epoch() << "Writing output" << endl;
	for(map<string, Scaf_stats>::const_iterator it = scaf_stats.begin(); it!=scaf_stats.end(); it++) {
		if(it->second.scaf_name != it->first) 
		{
			cerr << "Warning: " << it->first << " appears in mapping file(s) but not in assembly file" << endl;
		}
		else { 
			cout 	<< it->second.scaf_name << "\t" << it->second.size << "\t" << int(1000*it->second.gc)/10.0 << "\t" << int(100*it->second.coverage())/100.0 
				<< "\t" << it->second.Ns << endl;
		}
	}
	cerr << time_since_epoch() << "ok" << endl;
	cerr << time_since_epoch() << "Finished successfully" << endl << endl;
	
	return 0;
}

