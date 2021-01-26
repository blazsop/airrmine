BEGIN {
  RS="\n"
  FS="\t"
  OFS="\t"

  srand()

  # lc parameter can be specified from outside
  if( lc == "" ) linecount = 2500
  else           linecount = lc
}
BEGINFILE {
  # parse input repertoire filename and extract repertoire ID
  dbname = FILENAME
  sub(".*/",   "", dbname)     # remove directories from the full path
  sub(/\.tab/, "", dbname)     # remove ".tab" extension from end
  split( dbname, dbname_parts, "_" )
  repid = dbname_parts[2]

  total_lines = 0
  while( (getline < FILENAME) == 1 ) total_lines++
  close( FILENAME )

  if( linecount > total_lines )
	linecount = total_lines
  else
  {
    # generate "linecount" number of distinct random numbers
    for(i=1; i<=linecount; i++)
    {
      x = int(rand() * (total_lines-1)) + 2
      if( x in lines )
      {
        i--
        continue
      }
      else
        lines[ x ] = 1
    }
  }
}
FNR==1 { print $0 }
FNR>1  {
  if( linecount == total_lines || FNR in lines )
  {
    print $0
  }
}
ENDFILE {
  delete lines
}