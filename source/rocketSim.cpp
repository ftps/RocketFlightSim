#include "rocketSim.hpp"

RocketSim::RocketSim()
{

}

void RocketSim::timeMassThrustDragInput(const std::string& filename)
{
    std::fstream fs;
    fs.open(filename, std::fstream::in);

    if(!fs.is_open())
        std::cout<<"Unable to open file"<<"\n";


    std::vector<std::vector<double>> MainVec;
	std::vector<double> SubVec;
	std::string sFileRow,sFileCell;
	double dFileCell;
	size_t commaPos;
    

    while (std::getline(fs, sFileRow)) // extract line
	{
        if (sFileRow.find('#') != std::string::npos)
        {
            continue;
        }

		while (!sFileRow.empty())
		{
			commaPos = sFileRow.find(",");

			if (commaPos != std::string::npos)				// Process the Line
			{
				sFileCell = sFileRow.substr(0, commaPos);
				sFileRow.erase(0, commaPos + 1);
				dFileCell = std::stod(sFileCell);
				SubVec.push_back(dFileCell);
			}
			else
			{												// Process end of line ,pass to main vector, clear variables for loop
				dFileCell = std::stod(sFileRow);
				SubVec.push_back(dFileCell);
				MainVec.push_back(SubVec);
				SubVec.clear();
				sFileRow.clear();
			}
		}
	}

    for (size_t i = 0; i < MainVec.size(); i++)
    {   
        time.push_back(MainVec.at(i).at(0));
        mass.push_back(MainVec.at(i).at(5));
        thrust.push_back(MainVec.at(i).at(6));
        drag.push_back(MainVec.at(i).at(7));
    }
    

    fs.close();
}