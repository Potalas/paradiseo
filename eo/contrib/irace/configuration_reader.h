#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <string>
#include <memory>

#include <eo>
#include <ga.h>
#include <utils/checkpointing>
#include <eoInt.h>
#include <problems/eval/eoEvalIOH.h>
#include <filesystem>

#include <ioh.hpp>


using Ints = eoInt<eoMaximizingFitnessT<int>, int>;
using Bits = eoBit<eoMaximizingFitnessT<int>, int>;


/****************************************************************************
 configuration file format: csv using ";" as separator
 required columns: name, type, value, argument, probability
 exemple entries:
 
 instance_file;path
 
 ****************************************************************************/

/*
 A map of all files
 */

std::map<std::string, std::string> file_paths_library{
    {"instance_file", "instance.csv"},
    {"log_file", "onlymut_log.csv"}
};

/*==========================================================================
 ||                                                                       ||
 ||                                 Criterions                            ||
 ||                                                                       ||
 ==========================================================================*/

std::string criterion;

std::map<std::string, double> criteria{
    {"EAH", 0},
    {"NbGen", 0},
    {"Fitness", 0}
};

/*==========================================================================
 ||                                                                       ||
 ||                              parameters                               ||
 ||                                                                       ||
 ==========================================================================*/

/*
 Global sizes variables and a map with their names in the configuration file
 */


std::map<std::string, int> sizes_library {
    {"pop_size", 1},
    {"offspring_size", 2},
    {"nb_bucket_ECDF", 100},
    {"seed", int(time(0))},
    {"max_eval", 10000},
    {"max_generation", 1000},
    {"dimension", 0},
    {"instance",0},
    {"problem",0}
};

/*==========================================================================
 ||                                                                       ||
 ||                            Log switchs                                ||
 ||                                                                       ||
 ==========================================================================*/

/*
 Global booleans variables and a map with their names in the configuration file
 */

std::map<std::string, bool> switchs_library {
    {"gen_perf", true},   // if true then each generation performance is written in the logs
    {"best_perf", true}, // if true then at each generation the best found performance is written in the logs
    {"gen_sol", true},     // if true then each generation best solution is written in the logs
    {"best_sol", true},   // if true then at each generation the best found solution is written in the logs
    {"all_sol", false},     // if true then constructed solution is written in the logs
    {"output_mat", false},
    {"add_log", false},
    {"IOH_log", false}
};
/*==========================================================================
 ||                                                                       ||
 ||                        Crossover operators                            ||
 ||                                                                       ||
 ==========================================================================*/


eoForgeVector<eoQuadOp<Bits>> crossovers;
std::vector<double> crossovers_probabilities;


/*==========================================================================
 ||                                                                       ||
 ||                         Mutation operators                            ||
 ||                                                                       ||
 ==========================================================================*/

eoForgeVector<eoMonOp<Bits>> mutations;
std::vector<double> mutations_probabilities;

/*==========================================================================
 ||                                                                       ||
 ||                      Continuator operators                            ||
 ||                                                                       ||
 ==========================================================================*/

eoForgeVector<eoContinue<Bits>> continuators;
std::map<std::string, double> continuators_arguments
{
    {"FitContinue", 0},
    {"GenContinue", 0},
    {"EvalContinue", 0}
};

/*==========================================================================
 ||                                                                       ||
 ||                      Replacement operators                            ||
 ||                                                                       ||
 ==========================================================================*/

//Initialisation of all available replacements with their names in the configuration file
eoForgeVector<eoReplacement<Bits>> replacements;
std::vector<double> replacements_probabilities;


/*==========================================================================
 ||                                                                       ||
 ||                        Selector operators                             ||
 ||                                                                       ||
 ==========================================================================*/

eoForgeVector<eoSelectOne<Bits>> crossover_selectors;
eoForgeVector<eoSelectOne<Bits>> aftercross_selectors;
eoForgeVector<eoSelectOne<Bits>> mutation_selectors;

std::map<std::string, eoForgeVector<eoSelectOne<Bits>>> selectors
{
    {"crossover", crossover_selectors},
    {"aftercross", crossover_selectors},
    {"mutation", crossover_selectors}
};

std::vector<double> crossover_selectors_probabilities;
std::vector<double> aftercross_selectors_probabilities;
std::vector<double> mutation_selectors_probabilities;

std::map<std::string, std::vector<double>> selectors_probabilities
{
    {"crossover", crossover_selectors_probabilities},
    {"aftercross", aftercross_selectors_probabilities},
    {"mutation", mutation_selectors_probabilities}
};


/*-----------------------------------------------------------------------------------------------------------------*/

std::vector<std::string> split_line(std::string line, std::string sep = ";")
{
    std::vector<std::string> split_vector;
    size_t pos = 0;
    std::string token;
    while ((pos = line.find(sep)) != std::string::npos) {
        token = line.substr(0, pos);
        split_vector.push_back(token);
        line.erase(0, pos + sep.length());
    }
    split_vector.push_back(line);
    return split_vector;
}

std::vector<double> strings_to_doubles(std::vector<std::string> vect_string)
{
    std::vector<double> vect_double;
    for (std::string s : vect_string) {vect_double.push_back(std::stod(s));}
    return vect_double;
}

std::vector<int> strings_to_ints(std::vector<std::string> vect_string)
{
    std::vector<int> vect_int;
    for (std::string s : vect_string) {vect_int.push_back(std::stoi(s));}
    return vect_int;
}

std::vector<std::vector<int>> parse_sizes(std::string to_parse)
{
    std::vector<std::string> buckets_string = split_line(to_parse, ",");
    std::vector<std::vector<int>> buckets;
    for (std::string b : buckets_string)
    {
        buckets.push_back(strings_to_ints(split_line(b, ":")));
    }
    return buckets;
}

void read_configurations(std::string configuration_file_path)
{
    std::string sep = ";";
    std::ifstream configuration_file;
    configuration_file.open(std::filesystem::current_path().string() + "/" + configuration_file_path);
    /*
    if(configuration_file.is_open())
    {
        std::cout << std::filesystem::current_path().string() << std::endl;
        std::cout << "configuration file \"" << configuration_file_path << "\" doesn't exist" << std::endl;
    }
    else
    {
        std::cout << "configuration file \"" << configuration_file_path << "\" is open" << std::endl;
    }*/
    std::string line;
    std::map<std::string, int> columns{};
    
    // csv header
    configuration_file >> line;
    std::vector<std::string> row = split_line(line, sep);
    for (int column_index = 0; column_index < row.size(); column_index++) {columns[row[column_index]] = column_index;}
    
    std::string p_na, p_ty, p_def, p_val, p_arg, p_pro;
    
    std::vector<std::string> arguments;
    
    double cross_prob_total = 0;
    double mut_prob_total = 0;
    double remp_prob_total = 0;
    
    while(configuration_file >> line)
    {
        if (line.size() < 2)
        {
            break;
        }
        row = split_line(line, sep);
        
        p_na = row[columns["name"]];
        p_ty = row[columns["type"]];
        p_val = row[columns["value"]];
        p_arg = row[columns["argument"]];
        p_pro = row[columns["probability"]];
        
        if (p_ty == "path")
        {
            if (file_paths_library.find(p_na) != file_paths_library.end()) {file_paths_library[p_na] = p_val;}
            else {throw "ERROR: no " + p_na + " is required";}
        }
        if (p_ty == "parameter") {sizes_library[p_na] = std::stoi(p_val);}
        if (p_ty == "switch") {switchs_library[p_na] = (p_val == "true");}
        /* ---------------------------------- crossover ---------------------------------- */
        if (p_ty == "crossover")
        {
            // no argument
            if (p_na == "1PtBitXover")
            {crossovers.add<eo1PtBitXover<Bits>>();}
            
            // int argument
            if (p_na == "NPtsBitXover" )
            {crossovers.add<eoNPtsBitXover<Bits>>(std::stoi(p_arg));}
            
            // double argument
            if (p_na == "UBitXover" )
            {crossovers.add<eoUBitXover<Bits>>(std::stod(p_arg));}
            
            // int,int argument
            arguments = split_line(p_arg, ",");

            if (p_na == "BitGxOver")
            {crossovers.add<eoBitGxOver<Bits>>(std::stoi(arguments[0]), std::stoi(arguments[1]));}
            
            crossovers_probabilities.push_back(std::stod(p_pro));
            cross_prob_total += std::stod(p_pro);
            
        }
        /* ---------------------------------- mutation ---------------------------------- */
            
        if (p_ty == "mutation")
        {

            // no argument
            if (p_na == "OneBitFlip") {mutations.add<eoOneBitFlip<Bits>>();}
            if (p_na == "BitInversion") {mutations.add<eoBitInversion<Bits>>();}
            if (p_na == "BitNext") {mutations.add<eoBitNext<Bits>>();}
            if (p_na == "BitPrev") {mutations.add<eoBitPrev<Bits>>();}
            if (p_na == "UniformBitMutation") {mutations.add<eoUniformBitMutation<Bits>>();}
            
            // double argument
            if (p_na == "DetBitFlip") {mutations.add<eoDetBitFlip<Bits>>(std::stod(p_arg));}
            if (p_na == "DetSingleBitFlip") {mutations.add<eoDetSingleBitFlip<Bits>>(std::stod(p_arg));}
            if (p_na == "StandardBitMutation") {mutations.add<eoStandardBitMutation<Bits>>(std::stod(p_arg));}
            if (p_na == "ConditionalBitMutation") {mutations.add<eoConditionalBitMutation<Bits>>(std::stod(p_arg));}
            if (p_na == "ShiftedBitMutation") {mutations.add<eoShiftedBitMutation<Bits>>(std::stod(p_arg));}
            if (p_na == "FastBitMutation") {mutations.add<eoFastBitMutation<Bits>>(std::stod(p_arg));}
            
            // pair argument
            arguments = split_line(p_arg, ",");
            
            // double,bool argument
            if (p_na == "BitMutation") {mutations.add<eoBitMutation<Bits>>(std::stod(arguments[0]), arguments[1] == "true");}
            
            // double,double argument
            if (p_na == "NormalBitMutation") {mutations.add<eoNormalBitMutation<Bits>>(std::stod(arguments[0]), std::stod(arguments[1]));}
            
            // vector<vector<int>>|vector<double> argument
            arguments = split_line(p_arg, "|");
            if (p_na == "BucketBitMutation") {
                mutations.add<eoBucketBitMutation<Bits>>(parse_sizes(arguments[0]),strings_to_doubles(split_line(arguments[1], ",")));}
            
            mutations_probabilities.push_back(std::stod(p_pro));
            
            mut_prob_total += std::stod(p_pro);
        }

        /* ---------------------------------- continuator ---------------------------------- */
        if (p_ty == "continuator")
        {
            // int argument
            if (p_na == "GenContinue") {continuators.add<eoGenContinue<Bits>>(std::stoi(p_val));}
            //if (p_na == "EvalContinue") {continuators.add<eoEvalContinue<Bits>>(std::stoi(p_val));}
            
            // double argument
            if (p_na == "FitContinue") {continuators.add<eoFitContinue<Bits>>(std::stod(p_val));}
            
            continuators_arguments[p_na] = std::stod(p_val);
        }
        
        /* ---------------------------------- replacement ---------------------------------- */
        if (p_ty == "replacement")
        {

            // no argument
            if (p_na == "PlusReplacement") {replacements.add<eoPlusReplacement<Bits>>();}
            if (p_na == "CommaReplacement") {replacements.add<eoCommaReplacement<Bits>>();}
            if (p_na == "SSGAWorseReplacement") {replacements.add<eoSSGAWorseReplacement<Bits>>();}
            
            // int argument
            if (p_na == "SSGADetTournamentReplacement") {replacements.add<eoSSGADetTournamentReplacement<Bits>>(std::stoi(p_val));}
            if (p_na == "EPReplacement") {replacements.add<eoEPReplacement<Bits>>(std::stoi(p_val));}
            
            // double argument
            if (p_na == "SSGAStochTournamentReplacement") {replacements.add<eoSSGADetTournamentReplacement<Bits>>(std::stod(p_val));}

            replacements_probabilities.push_back(std::stod(p_pro));
            remp_prob_total += std::stod(p_pro);
        }
        
        /* ---------------------------------- selector ---------------------------------- */
        if (p_ty == "selector")
        {

            // no argument
            if (p_na == "RandomSelect") {selectors.at(p_val).add<eoRandomSelect<Bits>>();}
            if (p_na == "SequentialSelect") {selectors.at(p_val).add<eoSequentialSelect<Bits>>();}
            if (p_na == "ProportionalSelect") {selectors.at(p_val).add<eoProportionalSelect<Bits>>();}
            
            // int argument
            if (p_na == "DetTournamentSelect") {selectors.at(p_val).add<eoDetTournamentSelect<Bits>>(std::stoi(p_arg));}
            
            // double argument
            if (p_na == "StochTournamentSelect") {selectors.at(p_val).add<eoStochTournamentSelect<Bits>>(std::stod(p_arg));}
            
            selectors_probabilities.at(p_val).push_back(std::stod(p_pro));
        }
        
        /* ---------------------------------- criterion ---------------------------------- */
        if (p_ty == "criterion") {criterion = p_na;}
        
        configuration_file.ignore(1);

    }
    
    if (continuators.size() < 1) {throw "ERROR: required at least one continuator";}
    
    if (cross_prob_total > 1) {throw "ERROR: sum of crossovers probabilities (" + std::to_string(cross_prob_total) + ") must be <=  1";}
    crossovers_probabilities.push_back(1-cross_prob_total);
    
    if (mut_prob_total > 1) {throw "ERROR: sum of mutations probabilities (" + std::to_string(mut_prob_total) + ") must be <=  1";}
    mutations_probabilities.push_back(1-mut_prob_total);
    
    if (remp_prob_total > 1) {throw "ERROR: sum of remplacements probabilities (" + std::to_string(remp_prob_total) + ") must be <=  1";}
    
    for (auto sel : selectors)
    {
        if (sel.second.size() < 1)
        {
            selectors.at(sel.first).add<eoRandomSelect<Bits>>();
            selectors_probabilities.at(sel.first).push_back(1);
        }
    }
}

