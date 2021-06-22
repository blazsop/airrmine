# ------------------------------------------------------------------------------
# ImmChainTracer component - multifilter
#  - multiple germline-related sequence extractor  (v1.0 beta)
#
# input: standard AIRR database (.tsv)
# parameters:
#  - glfile: germline database file can be explicitly specified (optional)
#  - outdir: output directory
#  - jul_match: match JUNCTION_LENGTH (1 or 0)
#  - np1l_match: match NP1_LENGTH (1 or 0)
#  - np2l_match: match NP2_LENGTH (1 or 0)
# output: standard AIRR database files (.tsv) - current (or output) directory
#
#  by Peter Blazso
#  2021-06-21
# ------------------------------------------------------------------------------

BEGIN {
  # AIRR datafiles are tab-delimited BUT  (!)
  # because of a bug (?) in IMGT annotation or MakeDB.py database conversion
  # sometimes space(s) is/are left at the end of certain fields (like J_CALL)
  # this workaround applied
  FS = "\t"
  RS = "\n"

  # pre-defined column labels
  v_call_label = "v_call"
  j_call_label = "j_call"
  func_label   = "productive"
  seqid_label  = "sequence_id"
  jul_label    = "junction_length"
  np1l_label   = "np1_length"
  np2l_label   = "np2_length"

  # get database filename - if not given - exit with error
  dbfile = ARGV[1]
  if( dbfile == "" )
  {
    print "Please, specify a database (.tsv) file!" > "/dev/stderr"
    exit 1
  }

  # get or generate the filename containing relevant germline sequences
  if( glfile == "" )
  {
    split( dbfile, dbfile_parts, "_" )
    glfile = (dbfile_parts[1] "_germlines_" dbfile_parts[3] "_prepared.tsv")
  }

  # input parameter availability check and default settings
  # output directory (default: current "." directory)
  if( outdir == "" )   outdir = "."
  # jul_match (default: 1)
  if( jul_match == "" )  jul_match = 1
  # np1l_match (default: 0)
  if( np1l_match == "" ) np1l_match = 0
  # np2l_match (default: 0)
  if( np2l_match == "" ) np2l_match = 0

  # get and parse header line of germlines file store and copy
  getline < glfile
  header = $0


  # identify the indices of data columns
  for(i=1; i<=NF; i++)
  {
    if( $i == v_call_label ) v_call_idx = i
    if( $i == j_call_label ) j_call_idx = i
    if( $i == func_label   ) func_idx   = i
    if( $i == seqid_label  ) seqid_idx  = i
    if( $i == jul_label    ) jul_idx  = i
    if( $i == np1l_label   ) np1l_idx  = i
    if( $i == np2l_label   ) np2l_idx  = i
  }

  # build unique V_CALLs / J_CALLs arrays of FUNCTIONAL germline sequences
  # where index is the V_CALL or J_CALL and
  #       value is a comma-separated list of related SEQUENCE_IDs
  while( (getline < glfile) == 1 )
  {
    if( $func_idx == "TRUE" )
    {  
      # read, split v_call field and append v_calls list for each GL sequence
      split( $v_call_idx, v_call, "," )
      for(idx in v_call)
      {
        if( v_call[idx] in v_calls )
          v_calls[ v_call[idx] ] = v_calls[ v_call[idx] ] "," $seqid_idx
        else
          v_calls[ v_call[idx] ] = $seqid_idx
      }

      # read, split J_CALL field and append j_calls list for each GL sequence
      split( $j_call_idx, j_call, "," )
      for(idx in j_call)
      {
        if( j_call[idx] in j_calls )
          j_calls[ j_call[idx] ] = j_calls[ j_call[idx] ] "," $seqid_idx
        else
          j_calls[ j_call[idx] ] = $seqid_idx
      }

      # add the next SEQUENCE_ID to a list as an index
      seqids[ $seqid_idx ] = 0
      glseqs[ $seqid_idx ] = $0
    }
  }

  close( glfile )
}

BEGINFILE {
  # generate GL-related filelist based on input repertoire filename
  dbname = FILENAME
  sub(/\.tsv/, "", dbname)     # remove ".tsv" extension from end
  sub(".*/",   "", dbname)     # remove complete path from the beginning

  # add a header and the germline sequence
  for(idx in seqids)
  {
    outname = dbname "_" idx ".tsv"
    print header "\n" glseqs[ idx ] > (outdir "/" outname)
    print outname
  }
}

# match headers and if they are different, exit with error
FNR == 1 {
 if( $0 != header )
 {
   print "Germline and repertoire data headers mismatch!" > "/dev/stderr"
   exit 1
 }
}

# read and process input database file skipping header
FNR > 1 {
  # if the current sequence is non-functional just drop it
  if( $func_idx != "TRUE" ) next

  # read, split V_CALL field of db into "v_call" array and iterate through items
  split( $v_call_idx, v_call, "," )
  for( vi in v_call )
    if( v_call[vi] in v_calls )      # if current v_call has a match in GL list
    {                                # process J_CALLS the same way
      split( $j_call_idx, j_call, "," )
      for( ji in j_call )
        if( j_call[ji] in j_calls )  # if current j_call also has a match
        {
          # get and split the list of germlines belonging to the matched
          # v_call and j_call
          split( v_calls[v_call[vi]], v_gls, "," )
          split( j_calls[j_call[ji]], j_gls, "," )
          
          # mark all those germlines ONCE (in seqid array)
          # that are present in BOTH arrays  _AND_
          # share the same JUNCTION_LENGTH, NP1_LENGTH, NP2_LENGTH
          for( glvi in v_gls )
          for( glji in j_gls )          
            if( v_gls[glvi] == j_gls[glji] )
            {
              # get germline fields
              split( glseqs[ v_gls[glvi] ], gl_fields, "\t" )
              
              # mark ONLY if - in accordance with the "match" switches -
              # they also share the same JUNCTION_LENGTH, NP1, NP2
              if( $jul_idx  == gl_fields[ jul_idx ]  || !jul_match )
              if( $np1l_idx == gl_fields[ np1l_idx ] || !np1l_match )
              if( $np2l_idx == gl_fields[ np2l_idx ] || !np2l_match )          
                seqid[ v_gls[glvi] ] = 1
            }
        }
    }

  # write out the current sequence into all the marked, corresponding
  # germline-related target file
  for( gl in seqids )
    if( seqid[gl] == 1 )
    {
      seqid[gl] = 0
      print $0 >> (outdir "/" dbname "_" gl ".tsv")
    }
}
