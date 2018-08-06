// Logger

// write log to  file
void write_log_to_file(std::string _filename, Vcluster& _v_cl,double time, std::map<int,double> msrmt)
{
	// create_file for every processor involved
	std::ofstream output_file;
	std::string cpu_filename = _filename + "_" + std::to_string(_v_cl.getProcessUnitID());

	// check if the file exists:
	std::ifstream infile(cpu_filename);
	bool is_file_exists = infile.good();

	//open new one
	if (!is_file_exists)
	{
		//new file
		output_file.open (cpu_filename);
	}
	else
	{
		//Append
		output_file.open (cpu_filename, std::ios::app);
	}

	std::map<int,double>::iterator it_map;

	// if file doesn't exists write the header
	if(!is_file_exists)
	{
		output_file << "time" << "\t" ;

		for(it_map = msrmt.begin(); it_map != msrmt.end(); it_map++ )
		{
			output_file << it_map->first << "\t";
		}

		output_file << "\n";

	}

	// Write the collisions that occured this time step
	output_file << time << "\t" ;

	for(it_map = msrmt.begin(); it_map != msrmt.end(); it_map++ )
	{
		output_file << it_map->second << "\t";
		// Rest the collision counter
		it_map->second = 0;
	}

	output_file << "\n";

	output_file.close();


	return;
}
