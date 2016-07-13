# CPGobbler
A Python tool for extracting critical point and atomic basin information from ADF Tape21 files.

## To use
You'll need 1+ Tape21 files from ADF (not BAND), and an input text file we'll call input.txt.
You can use the tool on a machine that has Python 2.7+ with NumPy installed by executing the following command:

	python /path/to/CPGobbler.py /path/to/input.txt
 
The contents of input.txt will look something like this:

	/ADFData/Results.t21	[1,2,3,6,7,8]
	OtherResults.t21	[-1]
	*YetMoreResults*.t21	[1,2,3]

* Each line in input.txt contains a path to a Tape21 file and a list of critical point numbers to fetch.
  * The list of critical point numbers is comma-delimited without spaces, enclosed in square brackets `[ ]`.
  * The path to the Tape21 file and the list of critical point numbers are separated by a tab (that's **one tab**).
* The path to the file can be absolute as in line 1, or relative as in lines 2-3, and can also include wildcard (glob) characters.
  * If a relative path is provided then CPGobbler will assume the Tape21 file is in the same directory as input.txt, so on line 2, the file is assumed to be `/path/to/OtherResults.t21`.
  * Line 3 uses wildcard characters, where asterisk (\*) represents zero or more non-white space characters, and question mark (?) represents a single non-white space character.
  * Line 3 would return files such as
    * SP.LDA.YetMoreResults.t21 
    * 1\_YetMoreResults.t21 
    * YetMoreResults_GeometryOpt.t21 
    * 2016\_07\_13\_YetMoreResults.B3LYP.TZ2P.AllE.t21 

## Script options
At the top of CPGobbler.py, there are several blocks of data indicating which information should be collected for critical points and atomic basins:

	{'name':['x','y','z'],
	   'kf_name':'CP coordinates',
	   'collect':True},
	  
	  {'name':['rho'],
	   'kf_name':'CP density at',
	   'collect':True},
	  
	  {'name':['grad ' + i for i in ['x','y','z']],
	   'kf_name':'CP density gradient at',
	   'collect':True}
	   
### Enable/disable
The simplest option is to collect or not collect the data, which can be changed by setting the `True` to `False` or vice-versa.

### Adding new properties
There are three pieces of information necessary for each property:

#### `'name'`
A comma-delimited, bracketed list of strings for the name(s) of the data.
These lists can be populated manually, as with `['x','y','z']`, or through some other means such as a Python list comprehension ( `['grad ' + i for i in ['x','y','z']]` results in `['grad x','grad y',grad z']`, but it makes you feel cool to do it the complicated way!)
The composition of the name list is determined by the organization of the data for that particular property in the Tape21 file. For example, the XYZ coordinates are all contained in a single variable in the Tape21, as below

	x1 x2 x3
	x4 x5 x6
	y1 y2 y3
	y4 y5 y6
	z1 z2 z3
	z4 z5 z6
	
The number of items in the name list (3 for x,y,z) tells CPGobbler the separation of the data.
All the critical point data in a Tape21 is organized in this fashion, and you'll notice the same layout for properties to be fetched for atomic basins. 

If you need help with this, let me know.

#### `'kf_name'`
A string corresponding to the variable name in the Tape21 file.
Most of the informaiton in an output file is also in the Tape21.
To see what a given variable's kf name is, open the Tape21 file in ADFjobs by selecting the .t21 file in the job results and selecting SCM --> kf browser (SCM is now replaced by the SCM logo, and is the left-most button in the menu of ADF windows)

####`'collect'`
A boolean to determine whether a property should be collected.