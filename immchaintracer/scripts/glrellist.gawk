# ------------------------------------------------------------------------------
# ImmChainTracer component - glrellist
#  - germline-related file list generator  (v0.1 beta)
#
# input: AIRR database (NOT germline) file(s) (.tsv)
# output: list of germline-related AIRR database files
#
#  by Peter Blazso
# ------------------------------------------------------------------------------

BEGIN {
  # AIRR datafiles are tab-delimited BUT  (!)
  # because of a bug (?) in IMGT annotation or MakeDB.py database conversion
  # sometimes space(s) is/are left at the end of certain fields (like J_CALL)
  # this workaround applied
  FS = "\t"
  RS = "\n"

  # pre-defined column labels
  func_label   = "productive"
  seqid_label  = "sequence_id"

  # get database filename - if not given - exit with error
  dbfile = ARGV[1]
  if( dbfile == "" )
  {
    print "Please, specify a database (.tsv) file!" > "/dev/stderr"
    exit 1
  }
}

BEGINFILE {
  # generate GL-related filelist based on input repertoire filename
  dbdirs = FILENAME
  dbname = FILENAME
  sub(".*/",   "", dbname)     # remove directories from the full path
  sub(dbname,  "", dbdirs)     # remove database name from full path
  sub(/\.tsv/, "", dbname)     # remove ".tsv" extension from end

  split( dbname, dbname_parts, "_" )
  glfile = (dbdirs dbname_parts[1] "_germlines_" dbname_parts[3] \
            "_prepared.tsv")

  if( dbname_parts[2] == "germlines" ) nextfile

  # get and parse header line of germlines file store and copy
  getline < glfile
  header = $0

  # identify which columns contain FUNCTIONAL or SEQUENCE_ID strings
  for(i=1; i<=NF; i++)
  {
    if( $i == func_label   ) func_idx   = i
    if( $i == seqid_label  ) seqid_idx  = i
  }

  # print the list of germline-related AIRR database files
  # for the current repertoire database
  while( (getline < glfile) == 1 )
    if( $func_idx == "TRUE" )
    {  
      # add the next sequence_id to a list as an index
      printf( "%s ", (dbname "_" $seqid_idx ".tsv") )
    }

  close( glfile )
}
{
  nextfile
}