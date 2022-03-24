#include "configuration_reader.h"

class WModelFlat : public ioh::problem::wmodel::WModelOneMax
{
    public:
        WModelFlat(const int instance, const int n_variables,
                   const double dummy_para, const int epistasis_para, const int neutrality_para,
                   const int ruggedness_para)
        : WModelOneMax(instance, n_variables, dummy_para, epistasis_para, neutrality_para, ruggedness_para)
    {dimension = n_variables;}
        int dimension;
    protected:
        double transform_objectives(const double y) override
        { // Disable objective function shift & scaling.
            return y;
        }
};

WModelFlat read_instance(int problem_ID, int instance_ID, int instance_dimension)
{
    std::string sep = ";";
    std::string line;
    std::map<std::string, int> columns{};
    
    // csv header
    std::ifstream instance_file(std::filesystem::current_path().string() + "/" + file_paths_library["instance_file"]);
    
    instance_file >> line;
    std::vector<std::string> row = split_line(line, sep);
    for (int column_index = 0; column_index < row.size(); column_index++) {columns[row[column_index]] = column_index;}
    
    // instances characteristics
    do{
        instance_file >> line;
        row = split_line(line, sep);
    } while (std::stoi(row[columns["problem"]]) != problem_ID || std::stoi(row[columns["instance"]]) != instance_ID);

    if (instance_dimension == 0) {instance_dimension = std::stoi(row[columns["dimension"]]);}
    /*
    std::cout << "instance: " << std::stoi(row[columns["instance"]]) << std::endl;
    std::cout << "dimension: " << instance_dimension << std::endl;
    std::cout << "dummy: " << std::stod(row[columns["dummy"]]) << std::endl;
    std::cout << "epistasis: " << std::stoi(row[columns["epistasis"]]) << std::endl;
    std::cout << "neutrality: " << std::stoi(row[columns["neutrality"]]) << std::endl;
    std::cout << "ruggedness: " << std::stoi(row[columns["ruggedness"]]) << std::endl;
     */
    
    WModelFlat w_model_om(instance_ID,
                          instance_dimension,
                          std::stod(row[columns["dummy"]]),
                          std::stoi(row[columns["epistasis"]]),
                          std::stoi(row[columns["neutrality"]]),
                          std::stoi(row[columns["ruggedness"]]));
    
    return w_model_om;
}

/*
if (switchs_library["IOH_log"])
{
    ioh::trigger::OnImprovement on_improvement;
    ioh::watch::Evaluations evaluations;
    ioh::watch::TransformedYBest transformed_y_best;
    std::vector<std::reference_wrapper<ioh::logger::Trigger >> t = {on_improvement};
    std::vector<std::reference_wrapper<ioh::logger::Property>> w = {evaluations,transformed_y_best};
    csv_logger = std::make_shared<ioh::logger::FlatFile>(
        // {std::ref(on_improvement)},
        // {std::ref(evaluations),std::ref(transformed_y_best)},
        t, w,
        "GA",
        std::filesystem::create_directories("logs")
    );
    loggers.append(*csv_logger);
}
*/

