#include "instance_reader.h"


int main(int argc, char* argv[])
{
    read_configurations(argv[1]);
    WModelFlat w_model_om = read_instance(sizes_library["problem"], sizes_library["instance"], sizes_library["dimension"]);
    rng.reseed(sizes_library["seed"]);
    /*****************************************************************************
     * IOH stuff.
     *****************************************************************************/

    /***** IOH logger *****/
    if (continuators_arguments["EvalContinue"] == 0) {continuators_arguments["EvalContinue"] = continuators_arguments["GenContinue"] * (sizes_library["offspring_size"] + sizes_library["pop_size"]);}
    
    ioh::logger::eah::Log10Scale<double> target_range(0, continuators_arguments["FitContinue"], sizes_library["nb_bucket_ECDF"]);
    ioh::logger::eah::Log10Scale<size_t> budget_range(0, continuators_arguments["EvalContinue"], sizes_library["nb_bucket_ECDF"]);
    ioh::logger::EAH eah_logger(target_range, budget_range);

    ioh::logger::Combine loggers(eah_logger);

    std::shared_ptr<ioh::logger::FlatFile> csv_logger = nullptr;
    
    w_model_om.attach_logger(loggers);


    eoEvalIOHproblem<Bits> onemax_pb(w_model_om);
    
    eoBooleanGenerator<int> bgen;
    eoInitFixedLength<Bits> onemax_init(w_model_om.dimension, bgen);

    eoEvalCounterThrowException<Bits> eval_count(onemax_pb, continuators_arguments["EvalContinue"]);
    eoPopLoopEval<Bits> pop_eval(eval_count);
    
    eoPop<Bits> pop;
    pop.append(sizes_library["pop_size"], onemax_init);

    pop_eval(pop,pop);
    
    criteria["Fitness"] = pop.best_element().fitness();
    
    assert(pop.size() > 0);
    for(auto sol : pop) {assert(not sol.invalid());}

    // Set lambda to the pop size
    // if it was not set up at construction.
    if (sizes_library["offsprings_size"] == 0) {sizes_library["offsprings_size"] = pop.size();}

    
    std::ofstream logfile;
    logfile.open(file_paths_library["log_file"]);
    logfile << "generation;bestOffspring;bestOffspringValue;best;bestValue" << std::endl;

    criteria["NbGen"] = 0;

    int crossover_index;
    int mutation_index;
    int selector_index;
    int replacement_index;

    /*
    auto combconts = std::make_shared< std::vector<eoContinue<Bits>*> >();
    for (int c_index = 0 ; c_index < continuators.size(); c_index++)
    {
        
        combconts->push_back(&continuators.instantiate(c_index));
    }
    std::cout << combconts->size() << std::endl;*/
    if (w_model_om.objective().y < continuators_arguments["FitContinue"])
    {
        continuators_arguments["FitContinue"] = w_model_om.objective().y;
    }
    
    eoFunctorStore store;
    
    auto& fitcont = store.pack< eoFitContinue<Bits> >(continuators_arguments["FitContinue"]);
    auto& gencont = store.pack< eoGenContinue<Bits> >(continuators_arguments["GenContinue"]);
    auto combconts = std::make_shared< std::vector<eoContinue<Bits>*> >();
    combconts->push_back( &gencont );
    combconts->push_back( &fitcont );
    eoCombinedContinue<Bits> continuator(*combconts);
    
    // do while continuator is ok
    
    assert(pop.size() > 0);
    
    do {
        criteria["NbGen"]++;
        //std::cout << "seed " << sizes_library["seed"] << ": start generation " << criteria["NbGen"] << std::endl;

        eoPop<Bits> offsprings;
        while(offsprings.size() < sizes_library["offspring_size"]) {

            crossover_index = eo::rng.roulette_wheel(crossovers_probabilities);
            mutation_index = eo::rng.roulette_wheel(mutations_probabilities);

            if(crossover_index < crossovers.size()) {

                selector_index = eo::rng.roulette_wheel(selectors_probabilities["crossover"]);
                eoSelectOne<Bits>& crossover_selector = selectors["crossover"].instantiate(selector_index);
                crossover_selector.setup(pop);

                // Copy of const ref solutions,
                // because one alter them hereafter.
                Bits sol1 = crossover_selector(pop);
                Bits sol2 = crossover_selector(pop);

                // If the operator returns true,
                // solutions have been altered.
                eoQuadOp<Bits>& crossover = crossovers.instantiate(crossover_index);
                if(crossover(sol1, sol2)) {
                    sol1.invalidate();
                    sol2.invalidate();
                }

                // Select one of the two solutions
                // which have been crossed.
                eoPop<Bits> crossed;
                crossed.push_back(sol1);
                crossed.push_back(sol2);

                // The aftercross selector may need fitness,
                // so we evaluate those two solutions, if needed.
                pop_eval(crossed,crossed);
                
                selector_index = eo::rng.roulette_wheel(selectors_probabilities["aftercross"]);
                eoSelectOne<Bits>& aftercross_selector = selectors["aftercross"].instantiate(selector_index);
                aftercross_selector.setup(crossed);
                Bits sol3 = aftercross_selector(crossed);

                if(mutation_index < mutations.size())
                {
                    eoMonOp<Bits>& aftercross_mutation = mutations.instantiate(mutation_index);
                    if(aftercross_mutation(sol3)) {sol3.invalidate();}
                }
                offsprings.push_back(sol3);

            } else if (mutation_index < mutations.size())
            {
                selector_index = eo::rng.roulette_wheel(selectors_probabilities["mutation"]);
                eoSelectOne<Bits>& mutation_selector = selectors["mutation"].instantiate(selector_index);
                mutation_selector.setup(pop);
                Bits sol3 = mutation_selector(pop);
                
                eoMonOp<Bits>& mutation = mutations.instantiate(mutation_index);
                if(mutation(sol3)) {sol3.invalidate();}
                offsprings.push_back(sol3);
            }
        }
        assert(offsprings.size() == sizes_library["offsprings_size"]);
        //std::cout << "before eval " << continuators_arguments["EvalContinue"] << std::endl;
        pop_eval(pop, offsprings);
        //std::cout << "after eval" << std::endl;
        
        if (genStep == 1) {criteria["Fitness"] = pop.best_element().fitness();}
        if (criteria["Fitness"] < offsprings.best_element().fitness()) {criteria["Fitness"] = offsprings.best_element().fitness();}
        logfile << criteria["NbGen"] << ";";
        offsprings.best_element().printINlog(logfile);
        logfile << ";" << offsprings.best_element().fitness() << ";";
        
        replacement_index = eo::rng.roulette_wheel(replacements_probabilities);
        eoReplacement<Bits>& replacement = replacements.instantiate(replacement_index);
        //std::cout << "before replacement" << std::endl;
        replacement(pop, offsprings);
        
        pop.best_element().printINlog(logfile);
        logfile << ";" << criteria["Fitness"] << std::endl;
        //std::cout << "seed " << sizes_library["seed"] << ": end generation " << criteria["NbGen"] << std::endl;
    } while (continuator(pop));

    logfile.close();
    //assert(pop.size() > 0);
    for (auto sol : pop) {assert(not sol.invalid());}
               
    criteria["EAH"] = - ioh::logger::eah::stat::under_curve::volume(eah_logger);
    criteria["Fitness"] = - criteria["Fitness"];

    //std::cout << "seed " << sizes_library["seed"] << ": criterion " << criterion << std::endl;

    std::cout << criteria[criterion] << std::endl;
}
