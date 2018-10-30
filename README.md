# hmmscan_parser

Parser for hmmscan output from hmmer_3.2.1

Filters --domtblout output from hmmscan for results with an e-value of <1.0e-3 and coverage of >0.3.

Writes the domain name, domain length, query id, query length, e-value, hit start, hit end, query start, query end, and hit coverage to a tabular output file.
