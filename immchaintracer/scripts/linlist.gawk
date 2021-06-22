# ------------------------------------------------------------------------------
# ImmChainTracer component - linlist
#  - germline-related lineage file list generator  (v0.5 beta)
#
# input: AIRR germline file(s) (.tsv)
# output: list of germline-related lineage filenames
#
#  by Peter Blazso
# ------------------------------------------------------------------------------

BEGIN {
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
  dbname = FILENAME
  dbdirs = FILENAME
  sub(".*/",   "", dbname)     # remove directories from the full path
  sub(dbname,  "", dbdirs)     # remove database name from full path
  sub(/\.tsv/, "", dbname)     # remove ".tab" extension from end
  split( dbname, dbname_parts, "_" )

  # if not a germline-containing file is specified, jump over
  if( dbname_parts[2] != "germlines" ) nextfile

  # get and parse header line of germlines file store and copy
  getline < FILENAME
  header = $0

  # identify which columns contain FUNCTIONAL or SEQUENCE_ID strings
  for(i=1; i<=NF; i++)
  {
    if( $i == func_label   ) func_idx   = i
    if( $i == seqid_label  ) seqid_idx  = i
  }
  close( FILENAME )
}

FNR > 1 {
  # print the list of germline-related lineage names without suffix
  # for the current repertoire database
  if( $func_idx == "TRUE" )
  {  
    # add the next SEQUENCE_ID to a list as an index
    printf( "%s ", (dbname_parts[1] "_" $seqid_idx) )
  }
}
