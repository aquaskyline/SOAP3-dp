/*
 *
 *    IniParam.c
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#include "IniParam.h"


// This function is to load the input file for multi mode
MultiInputItem * loadMultiInputFile ( char * fileName, int isPair, int isReadBAM, int & lineNum )
{
    // return the array of MultiInputItem, and
    // the number of lines
    // get the number of '\n' inside the file
    lineNum = FUGetNumOfLines ( fileName );
    // create the array
    MultiInputItem * inputItemsArray = ( MultiInputItem * ) malloc ( ( lineNum + 1 ) * sizeof ( MultiInputItem ) );
    memset ( inputItemsArray, 0, ( lineNum + 1 ) * sizeof ( MultiInputItem ) );
    // read the file
    FILE * filein;
    char buffer[INPUT_BUFFER_SIZE];
    filein = ( FILE * ) fopen ( fileName, "r" );
    size_t bufferSize = fread ( buffer, sizeof ( char ), INPUT_BUFFER_SIZE, filein );
    int bufferIndex = 0;
    int currline = 0;
    int isEndOfLine = 0;
    char * tmpStr = ( char * ) malloc ( sizeof ( char ) * MAX_FIELD_LEN );

    while ( bufferSize > 0 && currline < lineNum )
    {
        // not end of file
        // get the query file 1
        FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                         inputItemsArray[currline].queryFile1,
                         isEndOfLine );

        if ( strlen ( inputItemsArray[currline].queryFile1 ) == 0 )
        {
            fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the name of query file 1.\n", fileName );
            exit ( 1 );
        }

        if ( isPair )
        {
            if ( isReadBAM == 0 )
            {
                // get the query file 2
                if ( !isEndOfLine )
                {
                    FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                     inputItemsArray[currline].queryFile2,
                                     isEndOfLine );
                }

                if ( strlen ( inputItemsArray[currline].queryFile2 ) == 0 )
                {
                    fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the name of query file 2.\n", fileName );
                    exit ( 1 );
                }
            }

            // get minimum value of insert size
            memset ( tmpStr, 0, sizeof ( char ) * MAX_FIELD_LEN ); // reset the array

            if ( !isEndOfLine )
            {
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 tmpStr, isEndOfLine );
            }

            if ( strlen ( tmpStr ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the minimum value of insert size.\n", fileName );
                exit ( 1 );
            }

            inputItemsArray[currline].insert_low = atoi ( tmpStr );
            // get maximum value of insert size
            memset ( tmpStr, 0, sizeof ( char ) * MAX_FIELD_LEN ); // reset the array

            if ( !isEndOfLine )
            {
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 tmpStr, isEndOfLine );
            }

            if ( strlen ( tmpStr ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the maximum value of insert size.\n", fileName );
                exit ( 1 );
            }

            inputItemsArray[currline].insert_high = atoi ( tmpStr );
        }

        // get the output prefix
        if ( !isEndOfLine )
        {
            FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                             inputItemsArray[currline].outputPrefix,
                             isEndOfLine );
        }

        if ( strlen ( inputItemsArray[currline].outputPrefix ) == 0 )
        {
            fprintf ( stderr, "Error. Inside %s, there exists a line which does not indicate the output prefix.\n", fileName );
            exit ( 1 );
        }

        if ( !isEndOfLine )
        {
            // get the readGrpID
            FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                             inputItemsArray[currline].readGrpID,
                             isEndOfLine );

            if ( strlen ( inputItemsArray[currline].readGrpID ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which contains an empty field for read group ID.\n", fileName );
                exit ( 1 );
            }

            // get the sampleName
            if ( !isEndOfLine )
            {
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 inputItemsArray[currline].sampleName,
                                 isEndOfLine );
            }

            if ( strlen ( inputItemsArray[currline].sampleName ) == 0 )
            {
                fprintf ( stderr, "Error. Inside %s, there exists a line which contains read group ID but does not indicate the sample name.\n", fileName );
                exit ( 1 );
            }

            updateUserInputText ( inputItemsArray[currline].sampleName );

            if ( !isEndOfLine )
            {
                // get the readGrpOpt
                FUGetNextField ( filein, buffer, bufferSize, bufferIndex,
                                 inputItemsArray[currline].readGrpOpt,
                                 isEndOfLine );

                if ( strlen ( inputItemsArray[currline].readGrpOpt ) == 0 )
                {
                    fprintf ( stderr, "Error. Inside %s, there exists a line which contains an empty field for read group other option.\n", fileName );
                    exit ( 1 );
                }

                updateUserInputText ( inputItemsArray[currline].readGrpOpt );
            }
        }

        if ( !isEndOfLine )
        {
            // skip all characters until the end of line
            FUSkipToEOL ( filein, buffer, bufferSize, bufferIndex );
        }

        currline++;
    }

    if ( currline == 0 )
    {
        fprintf ( stderr, "Error. The file %s is empty!", fileName );
        exit ( 1 );
    }

    fclose ( filein );
    return inputItemsArray;
}

// This function replace the special character \n or \t with a newline or a tab
void updateUserInputText ( char * str )
{
    // "\t" -> '\t'
    // "\n" => '\n'
    int i, j;

    for ( i = 0; i < strlen ( str ) - 1; i++ )
    {
        if ( str[i] == '\\' )
        {
            if ( str[i + 1] == 't' )
            {
                str[i] = '\t';

                // move all the remaining characters one-pace forward
                for ( j = i + 1; j < strlen ( str ) - 1; j++ )
                {
                    str[j] = str[j + 1];
                }

                str[strlen ( str ) - 1] = '\0';
            }
            else if ( str[i + 1] == 'n' )
            {
                str[i] = '\n';

                // move all the remaining characters one-pace forward
                for ( j = i + 1; j < strlen ( str ) - 1; j++ )
                {
                    str[j] = str[j + 1];
                }

                str[strlen ( str ) - 1] = '\0';
            }
        }
    }
}

// This function is to parse the ini file and collect the paramaters of
// 1. The file extension of SA file
// 2. The number of CPU threads
// 3. The alignment model
// 4. The memory size in the GPU card
int ParseIniFile ( char * iniFileName, IniParams & ini_params )
{
    dictionary * ini;
    char tmp[10];
    ini = iniparser_load ( iniFileName, FALSE );

    if ( ini == NULL )
    {
        printf ( "[ParseIniFile] File not found!\n" );
        return -1;
    }

    // Database:SaValueFileExt parameter
    iniparser_copystring ( ini, "Database:SaValueFileExt", ini_params.Ini_SaValueFileExt, ".sa", MAX_FILEEXT_LEN );
    // Alignment:NumOfCpuThreads parameter
    ini_params.Ini_NumOfCpuThreads = iniparser_getuint ( ini, "Alignment:NumOfCpuThreads", 3 );

    if ( ini_params.Ini_NumOfCpuThreads <= 0 || ini_params.Ini_NumOfCpuThreads > MAX_NUM_CPU_THREADS )
    {
        printf ( "[ParseIniFile] Invalid value for Ini_NumOfCpuThreads(%u)!\n", ini_params.Ini_NumOfCpuThreads );
        return -1;
    }

    // Alignment:HostAlignmentModel parameter
    iniparser_copystring ( ini, "Alignment:HostAlignmentModel", ini_params.Ini_HostAlignmentModelStr, "16G", 3 );

    if ( strcmp ( ini_params.Ini_HostAlignmentModelStr, "8G" ) == 0 )
    {
        ini_params.Ini_HostAlignmentModel = SRA_MODEL_8G;
    }
    else if ( strcmp ( ini_params.Ini_HostAlignmentModelStr, "16G" ) == 0 )
    {
        ini_params.Ini_HostAlignmentModel = SRA_MODEL_16G;
    }
    else
    {
        printf ( "[ParseIniFile] Invalid value for HostAlignmentModel(%s))!\n", ini_params.Ini_HostAlignmentModelStr );
        return -1;
    }

    // Alignment:Soap3MisMatchAllow parameter
    ini_params.Ini_Soap3MisMatchAllow = iniparser_getint ( ini, "Alignment:Soap3MisMatchAllow", 2 );

    if ( ini_params.Ini_Soap3MisMatchAllow < 0 || ini_params.Ini_Soap3MisMatchAllow > 4 )
    {
        printf ( "[ParseIniFile] The program only supports number of mismatches for soap3 alignment between 0 and 4.\n" );
        return -1;
    }

    // DP:DPScoreThreshold parameter
    iniparser_copystring ( ini, "DP:DPScoreThreshold", tmp, "DEFAULT", 10 );

    if ( strcmp ( tmp, "DEFAULT" ) == 0 )
    {
        ini_params.Ini_isDefaultThreshold = 1;
    }
    else
    {
        ini_params.Ini_isDefaultThreshold = 0;
        ini_params.Ini_DPScoreThreshold = atoi ( tmp );
    }

    ini_params.Ini_MaxOutputPerRead = iniparser_getuint ( ini, "Alignment:MaxOutputPerRead", 1000 );
    // Device:GPUMemory parameter
    ini_params.Ini_GPUMemory = iniparser_getint ( ini, "Device:GPUMemory", 3 );

    if ( ini_params.Ini_GPUMemory < 3 )
    {
        printf ( "[ParseIniFile] The program only supports GPU with at least 3G memory.\n" );
        return -1;
    }

    iniparser_copystring ( ini, "PairEnd:StrandArrangement", tmp, "+/-", 10 );

    if ( strcmp ( tmp, "+/+" ) == 0 )
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_POS_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_POS_STRAND;
    }
    else if ( strcmp ( tmp, "-/+" ) == 0 )
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_NEG_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_POS_STRAND;
    }
    else if ( strcmp ( tmp, "-/-" ) == 0 )
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_NEG_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_NEG_STRAND;
    }
    else
    {
        ini_params.Ini_PEStrandLeftLeg = QUERY_POS_STRAND;
        ini_params.Ini_PEStrandRightLeg = QUERY_NEG_STRAND;
    }

#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
    ini_params.Ini_PEStrandLeftLeg = QUERY_POS_STRAND;
    ini_params.Ini_PEStrandRightLeg = QUERY_POS_STRAND;
#endif
    ini_params.Ini_PEMaxOutputPerPair = iniparser_getuint ( ini, "PairEnd:MaxOutputPerPair", 1000 );
    ini_params.Ini_MaxHitsEachEndForPairing = iniparser_getuint ( ini, "PairEnd:MaxHitsEachEndForPairing", 10000 );
    // FOR DP MODULE
    ini_params.Ini_MatchScore = iniparser_getint ( ini, "DP:MatchScore", 1 );
    ini_params.Ini_MismatchScore = iniparser_getint ( ini, "DP:MismatchScore", -2 );
    ini_params.Ini_GapOpenScore = iniparser_getint ( ini, "DP:GapOpenScore", -3 );
    ini_params.Ini_GapExtendScore = iniparser_getint ( ini, "DP:GapExtendScore", -1 );
    // score
    ini_params.Ini_minMAPQ = iniparser_getint ( ini, "Score:MinMAPQ", 1 );
    ini_params.Ini_maxMAPQ = iniparser_getint ( ini, "Score:MaxMAPQ", 40 );
    // index to be shared among multiple copies of soap3-dp
    ini_params.Ini_shareIndex = iniparser_getint ( ini, "Index:ShareIndex", 0 );
    // max length allowed for read name
    ini_params.Ini_maxReadNameLen = iniparser_getint ( ini, "Reads:MaxLenReadName", 64 );
    // max length from the front of the read for clipping
    ini_params.Ini_maxFrontLenClipped = iniparser_getint ( ini, "Clipping:MaxFrontLenClipped", 3 );
    // max length from the end of the read for clipping
    ini_params.Ini_maxEndLenClipped = iniparser_getint ( ini, "Clipping:MaxEndLenClipped", 8 );
    // whether the seed will proceed to perform DP if there are too many hits
    ini_params.Ini_proceedDPForTooManyHits = iniparser_getint ( ini, "OtherSettings:ProceedDPForTooManyHits", 0 );
    // whether the read will perform SOAP3 module
    ini_params.Ini_skipSOAP3Alignment = iniparser_getint ( ini, "OtherSettings:SkipSOAP3Alignment", 0 );
    // whether the bwa-like MAPQ score should be reported
    ini_params.Ini_bwaLikeScore = iniparser_getint ( ini, "Score:BWALikeScore", 0 );
    iniparser_freedict ( ini );
    return 0;
}


void INIPrintAlignmentUsage ( InputOptions * inputOptions, char * programName )
{
    if ( inputOptions->readType == SINGLE_READ )
    {
        if ( !inputOptions->isReadList )
        {
            UIprintUsageSingle ( programName );
        }
        else
        {
            UIprintUsageSingleList ( programName );
        }
    }
    else if ( inputOptions->readType == PAIR_END_READ )
    {
        if ( !inputOptions->isReadList )
        {
            UIprintUsagePair ( programName );
        }
        else
        {
            UIprintUsagePairList ( programName );
        }
    }
    else
    {
        UIPrintUsageOverview ( programName );
    }
}

bool parseInputArgs ( int argc, char ** argv, InputOptions & input_options )
{
    char * programName = argv[0];
    bool readGroupSupplied = false; // whether an option "-D" is specified
    bool sampleNameSupplied = false; // whether an option "-A" is specified
    bool readGroupOptionSupplied = false; // whether an option "-R" is specified

    // This function is to check and parse the input arguments
    // arguments: <program> single bwtCodeIndex queryFileName numQueries maxReadLength [options]
    // OR         <program> pair bwtCodeIndex queryFileName1 queryFileName2 numQueries maxReadLength [options]
    // If the argument is not correct, then it returns FALSE. Else it returns TRUE
    if ( argc < 2 )
    {
        // fprintf(stderr,"No argument has been provided.\n");
        // printUsage(argv[0]);
        UIPrintUsageOverview ( programName );
        return false;
    }

    if ( strcmp ( argv[1], "single" ) == 0 )
    {
        input_options.readType = SINGLE_READ;
        input_options.isReadList = FALSE;
    }
    else if ( strcmp ( argv[1], "pair" ) == 0 )
    {
        input_options.readType = PAIR_END_READ;
        input_options.isReadList = FALSE;
    }
    else if ( strcmp ( argv[1], "single-multi" ) == 0 )
    {
        input_options.readType = SINGLE_READ;
        input_options.isReadList = TRUE;
    }
    else if ( strcmp ( argv[1], "pair-multi" ) == 0 )
    {
        input_options.readType = PAIR_END_READ;
        input_options.isReadList = TRUE; // list of read files
    }
    else
    {
        fprintf ( stderr, "Miss the keyword 'single', 'single-multi', 'pair' or 'pair-multi'.\n" );
        return false;
    }

    if ( argc == 2 )
    {
        INIPrintAlignmentUsage ( &input_options, programName );
        return false;
    }

    // check whether the input read type is BAM
    input_options.isReadBAM = 0;

    for ( int i = 2; i < argc; i++ )
    {
        if ( strcmp ( argv[i], "-bam" ) == 0 || strcmp ( argv[i], "-BAM" ) == 0 )
        {
            input_options.isReadBAM = 1;
            break;
        }
    }

    // check the number of arguments for the specific read type
    int min_num_args = ( input_options.readType == PAIR_END_READ && input_options.isReadBAM == 0 && input_options.isReadList == 0 ) ? 5 : 4;

    if ( argc < min_num_args )
    {
        fprintf ( stderr, "Invalid number of command-line arguments.\n" );
        INIPrintAlignmentUsage ( &input_options, programName );
        return false;
    }

    input_options.indexName = argv[2];
    // get the query file name(s) and the maximum read length
    input_options.queryFileName = argv[3];

    // get the second query file name
    // only for "pair" mode and not a BAM input
    if ( input_options.readType == PAIR_END_READ && !input_options.isReadList && !input_options.isReadBAM )
    {
        input_options.queryFileName2 = argv[4];
    }

    // set default values
    input_options.maxReadLength = DEFAULT_MAX_READ_LEN;
    input_options.numMismatch = -1;
    input_options.insert_low = 1;
    input_options.insert_high = 500;
    input_options.outputFormat = SRA_OUTPUT_FORMAT_SAM_API;
    input_options.isOutputBinary = 0;
    input_options.outputPrefix = argv[3];
    input_options.enableDP = 1;
    input_options.isIlluminaQual = 0;
    input_options.GPUDeviceID = -1;
    input_options.readGroup = ( char * ) "";
    input_options.sampleName = ( char * ) "";
    input_options.readGrpOption = ( char * ) ""; // read group option
    // set default alignment type = all-best alignment
    input_options.alignmentType = OUTPUT_ALL_BEST;
    input_options.isPrintMDNM = false;

    // get the options
    for ( int i = min_num_args; i < argc; i++ )
    {
        if ( strcmp ( argv[i], "-h" ) == 0 )
        {
            // output option
            // 1: all valid alignments
            // 2: all best alignments (default)
            // 3: unique best alignments
            // 4: random best alignments
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the output option after '-h' ( 1, 2, 3 or 4 )\n" );
                return false;
            }

            int input_type = atoi ( argv[++i] );

            if ( input_type == 1 )
            {
                input_options.alignmentType = OUTPUT_ALL_VALID;
            }
            else if ( input_type == 2 )
            {
                input_options.alignmentType = OUTPUT_ALL_BEST;
            }
            else if ( input_type == 3 )
            {
                input_options.alignmentType = OUTPUT_UNIQUE_BEST;
            }
            else if ( input_type == 4 )
            {
                input_options.alignmentType = OUTPUT_RANDOM_BEST;
            }
            else
            {
                fprintf ( stderr, "The output option should be 1, 2, 3 or 4\n" );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-l" ) == 0 || strcmp ( argv[i], "-L" ) == 0 )
        {
            // the length of the longest read in the input
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the length after '-L'\n" );
                return false;
            }

            input_options.maxReadLength = atoi ( argv[++i] );

            if ( input_options.maxReadLength < 0 )
            {
                fprintf ( stderr, "The length should not be less than 0\n" );
                return false;
            }
            else if ( input_options.maxReadLength > MAX_READ_LENGTH )
            {
                fprintf ( stderr, "The length should not be greater than %u\n", MAX_READ_LENGTH );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-u" ) == 0 )
        {
            // the upper bound of the insert size
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the maximum value of insert size after '-u'\n" );
                return false;
            }

            input_options.insert_high = atoi ( argv[++i] );
        }
        else if ( strcmp ( argv[i], "-v" ) == 0 )
        {
            // the lower bound of the insert size
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the minimum value of insert size after '-v'\n" );
                return false;
            }

            input_options.insert_low = atoi ( argv[++i] );
        }
        else if ( strcmp ( argv[i], "-b" ) == 0 )
        {
            // the output format
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the output file format after '-b'\n" );
                return false;
            }

            input_options.outputFormat = atoi ( argv[++i] );

            if ( input_options.outputFormat != SRA_OUTPUT_FORMAT_PLAIN &&
                    input_options.outputFormat != SRA_OUTPUT_FORMAT_SAM_API &&
                    input_options.outputFormat != SRA_OUTPUT_FORMAT_BAM )
            {
                fprintf ( stderr, "Incorrect Output Format supplied. Please refer to Usage." );
                return false;
            }

            if ( input_options.outputFormat == SRA_OUTPUT_FORMAT_BAM )
            {
                input_options.outputFormat = SRA_OUTPUT_FORMAT_SAM_API;
                input_options.isOutputBinary = 1;
            }
        }
        else if ( strcmp ( argv[i], "-s" ) == 0 )
        {
            // disable dynamic programming for the unmapped mate
            input_options.enableDP = 0;

            if ( i + 1 < argc && isdigit ( argv[i + 1][0] ) )
            {
                int tmp = atoi ( argv[i + 1] );

                if ( tmp >= 0 && tmp <= 4 )
                {
                    input_options.numMismatch = tmp;
                    i++;
                }
                else
                {
                    fprintf ( stderr, "The value after '-s' should be between 0 and 4\n" );
                }
            }
        }
        else if ( strcmp ( argv[i], "-o" ) == 0 )
        {
            // specify the output file prefix
            if ( input_options.isReadList == 1 )
            {
                fprintf ( stderr, "For multiple sets of input files, the output file prefix could not be specified." );
                fprintf ( stderr, "The option '-o' is thus ignored." );
            }
            else
            {
                if ( i + 1 >= argc )
                {
                    fprintf ( stderr, "Please specify the output file prefix after '-o'\n" );
                    return false;
                }

                input_options.outputPrefix = argv[++i];
            }
        }
        else if ( strcmp ( argv[i], "-I" ) == 0 )
        {
            // Quality is in Phred+64 format
            input_options.isIlluminaQual = 1;
        }
        else if ( strcmp ( argv[i], "-D" ) == 0 )
        {
            // to assign the read group ID
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the read group ID after '-D'\n" );
                return false;
            }

            input_options.readGroup = argv[++i];
            readGroupSupplied = true;
        }
        else if ( strcmp ( argv[i], "-A" ) == 0 )
        {
            // to assign the sample name
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the sample name after '-A'\n" );
                return false;
            }

            input_options.sampleName = argv[++i];
            // update the input text
            updateUserInputText ( input_options.sampleName );
            sampleNameSupplied = true;
        }
        else if ( strcmp ( argv[i], "-R" ) == 0 )
        {
            // to assign the read group option
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the sample name after '-R'\n" );
                return false;
            }

            input_options.readGrpOption = argv[++i];
            // update the input text
            updateUserInputText ( input_options.readGrpOption );
            readGroupOptionSupplied = true;
        }
        else if ( strcmp ( argv[i], "-c" ) == 0 )
        {
            // to specific the GPU device ID
            if ( i + 1 >= argc )
            {
                fprintf ( stderr, "Please specify the GPU device ID after '-c'\n" );
                return false;
            }

            input_options.GPUDeviceID = atoi ( argv[++i] );

            if ( input_options.GPUDeviceID < 0 )
            {
                fprintf ( stderr, "The GPU device ID should not be less than 0\n" );
                return false;
            }
        }
        else if ( strcmp ( argv[i], "-p" ) == 0 )
        {
            // print MD string and NM tag
            input_options.isPrintMDNM = true;
        }
        else 
        {
            fprintf ( stderr, "Unknown option \"%s\". Please refer to Usage.\n", argv[i]);
            return false;
        }
    }

    if ( input_options.readType == PAIR_END_READ )
    {
        if ( input_options.insert_high == -1 )
        {
            fprintf ( stderr, "Please specify the maximum value of insert size.\n" );
            return false;
        }

        if ( input_options.insert_low == -1 )
        {
            fprintf ( stderr, "Please specify the minimum value of insert size.\n" );
            return false;
        }

        if ( input_options.insert_low > input_options.insert_high )
        {
            fprintf ( stderr, "The minimum value of insert size should not be greater than the maximum value of insert size.\n" );
            return false;
        }
    }

    // only allow three cases:
    // case 1: user does not specify "-D", "-A" or "-R"
    // case 2: user only specifies "-D" and "-A"
    // case 3: user only specifies "-D", "-A" and "-R"
    if ( readGroupOptionSupplied && ( ( !sampleNameSupplied ) || ( !readGroupSupplied ) ) )
    {
        fprintf ( stderr, "Since an option '-D' is specified, both options '-A' and '-R' have to be specified too.\n" );
        return false;
    }

    if ( sampleNameSupplied && ( !readGroupSupplied ) )
    {
        fprintf ( stderr, "Since an option '-A' is specified, the options '-D' has to be specified too.\n" );
        return false;
    }

    if ( readGroupSupplied && ( !sampleNameSupplied ) )
    {
        fprintf ( stderr, "Since an option '-D' is specified, the options '-A' has to be specified too.\n" );
        return false;
    }

    return true;
}


// print out the parameters
void printParameters ( InputOptions input_options, IniParams ini_params )
{
    printf ( "\n----------Parameters------------------\n" );
    printf ( "Number of mismatches: %i\n", input_options.numMismatch );

    if ( input_options.alignmentType == OUTPUT_ALL_VALID )
    { printf ( "Output all valid alignments\n" ); }
    else if ( input_options.alignmentType == OUTPUT_ALL_BEST )
    { printf ( "Ouput all best alignments\n" ); }
    else if ( input_options.alignmentType == OUTPUT_UNIQUE_BEST )
    { printf ( "Ouput unique best alignments\n" ); }
    else if ( input_options.alignmentType == OUTPUT_RANDOM_BEST )
    { printf ( "Ouput random best alignments\n" ); }

    printf ( "Number of CPU threads: %u\n", ini_params.Ini_NumOfCpuThreads );

    if ( input_options.insert_high > 0 )
    { printf ( "Upper bound of insert size: %u\n", input_options.insert_high ); }

    if ( input_options.insert_low > 0 )
    { printf ( "Lower bound of insert size: %u\n", input_options.insert_low ); }

    printf ( "enableDP: %i\n", ( int ) input_options.enableDP );
    printf ( "outputPrefix: %s\n", input_options.outputPrefix );
    printf ( "readGroup: %s\n", input_options.readGroup );
    printf ( "readGrpOption: %s\n", input_options.readGrpOption );
    printf ( "--------------------------------------\n\n" );
}

// print out the parameters
void printDPParameters ( DPParameters dp_params )
{
#ifdef BGS_OUTPUT_PARAMETERS
    printf ( "\n----------DP Parameters------------------\n" );
    printf ( "Match score: %i\n", dp_params.matchScore );
    printf ( "Mismatch score: %i\n", dp_params.mismatchScore );
    printf ( "Open gap score: %i\n", dp_params.openGapScore );
    printf ( "Extend gap score: %i\n", dp_params.extendGapScore );
    printf ( "numOfCPUThreads: %i\n", dp_params.numOfCPUThreads );
    printf ( "numOfCPUForSeeding: %i\n", dp_params.numOfCPUForSeeding );
    printf ( "softClipLeft: %i\n", dp_params.softClipLeft );
    printf ( "softClipRight: %i\n", dp_params.softClipRight );
    printf ( "tailTrimLen: %i\n", dp_params.tailTrimLen );
    printf ( "singleDPSeedNum: %i\n", dp_params.singleDPSeedNum );

    for ( int i = 0; i < 3; i++ )
    {
        printf ( "singleDPSeedPos[%i]: %i\n", i, dp_params.singleDPSeedPos[i] );
    }

    for ( int i = 0; i < 2; i++ )
    {
        printf ( "cuffoffthreshold[%i]: %i\n", i, dp_params.paramRead[i].cutoffThreshold );
        printf ( "maxHitNum[%i]: %i\n", i, dp_params.paramRead[i].maxHitNum );
        printf ( "sampleDist[%i]: %i\n", i, dp_params.paramRead[i].sampleDist );
        printf ( "seedLength[%i]: %i\n", i, dp_params.paramRead[i].seedLength );
    }

    printf ( "--------------------------------------\n\n" );
#endif
}


// This function is to update the values inside InputOptions
// according to the i-th set of values inside the array of MultiInputItem
void updateInputOption ( InputOptions * input_options, MultiInputItem * multiInputArray, int i )
{
    input_options->queryFileName = multiInputArray[i].queryFile1;
    input_options->queryFileName2 = multiInputArray[i].queryFile2;
    input_options->insert_low = multiInputArray[i].insert_low;
    input_options->insert_high = multiInputArray[i].insert_high;
    input_options->outputPrefix = multiInputArray[i].outputPrefix;
    input_options->readGroup = multiInputArray[i].readGrpID;
    input_options->sampleName = multiInputArray[i].sampleName;
    input_options->readGrpOption = multiInputArray[i].readGrpOpt;
}
