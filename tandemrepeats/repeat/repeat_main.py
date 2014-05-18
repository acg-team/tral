import sys, logging

from tandemrepeats.repeat import repeat, repeat_io

def worker(parameters = None):

    ''' SET UP YOUR PARAMETER OPTIONS '''

    # Specify which scores shall be calculated.
    # Options are: ['phylo', 'phylo_gap01', 'phylo_gap001', 'entropy', 'parsimony','pSim']
    scoring_functions = ['phylo', 'phylo_gap01', 'phylo_gap001']

    # Specify the path where the repeat is saved (<path/to/file>)
    # Specify the repeat's sequence_type ("AA" or "DNA")
    repeats = repeat_io.read_repeats('repeat_dna.txt', sequence_type="DNA")

    # If the output should be save to a file, define <output_file_path>.
    # Else, set output_file_path = None
    output_file_path = None
    if output_file_path == None:
        print_function = sys.stdout.write
    else:
        f = open(output_file_path, 'w')
        print_function = f.write

    scoring_functions_displayed = []
    for iC in scoring_functions:
        if iC in ['phylo', 'phylo_gap01', 'phylo_gap001']:
            scoring_functions_displayed.extend([iC, iC + "_ignore_coherent_deletions", iC + "_ignore_trailing_gaps", iC + "_ignore_trailing_gaps_and_coherent_deletions"])
        else:
            scoring_functions_displayed.append(iC)



    for iSeq,iR in repeats.items():

        iR.calculate_scores(scoreslist=scoring_functions)
        iR.calculate_pValues(scoreslist=scoring_functions)

        # Print the ID of the Repeat (Comment/Uncomment)
        print_function(iSeq + '\n')

        # Print scoring_function names (Comment/Uncomment)
        print_function("\t".join(scoring_functions_displayed)+ '\n')

        # Print the corresponding pValues (Comment/Uncomment)
        show_pValues = [str(iR.pValue[i]) for i in scoring_functions_displayed]
        print_function("\t".join(show_pValues)+ '\n')

        # Print the corresponding divergences (Comment/Uncomment)
        show_divergences = [str(iR.divergence[i]) for i in scoring_functions_displayed]
        print_function("\t".join(show_divergences)+ '\n\n')


    if not output_file_path == None:
        f.close()

if __name__=="__main__":

    #worker(jobID = sys.argv[1])
    worker()
