convertListToS4               <-  function( list , class ) {  
  S4Object                    <-  new( class ) 
  slotNames                   <-  slotNames( S4Object )
  listNames                   <-  names( list )
  checkNames                  <-  all( listNames   %in%  slotNames ) 
  checkEqualNameLength        <-  length( listNames ) == length( slotNames  )
  if( ! checkNames ) {
    stop( "names of 'list' do not match 'class'" )
  }
  if( ! checkEqualNameLength ) {
    warning( "some elements will be lost by conversion to S4 object" )
  }
  for( slot in slotNames) {
    slot( S4Object  , slot, check = FALSE ) <-  list[[ slot ]]  # check whole object at once 
  }
  isValidObject               <-  validObject( S4Object )
  return( S4Object )	    
}
