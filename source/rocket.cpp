#include "rocket.hpp"

Component::Component(const std::string& aux_s, const double& Xtip, const uint& type) : type(type)
{
	std::vector<double> aux_v;

	// DARK MAGIC TO GET VALUES FROM STRING
	for(const std::string& s : spltString(aux_s)){
		aux_v.emplace_back(std::stod(s));
	}

	// Component Geometry
	L = aux_v[0];
	R = aux_v[1];
	M = aux_v[2];
	Xcm = Xtip + ((type == cone) ? 2.0/3.0 : 1.0/2.0)*L;
}

fm::Matrix33 Component::momInertia()
{
	double Ixx, Izz;

	// Cone
	if(type == cone){
		Ixx = (1.0/10.0) * M * fm::sq(L) + (3.0/20.0) * M * fm::sq(R);
    	Izz = (3.0/10.0) * M * fm::sq(R);
	} // Cylinder
	else{
		Ixx = (1.0/12.0) * M * (3 * fm::sq(R) + fm::sq(L));
		Izz = (1.0/2.0) * M * fm::sq(R);
	}

	return fm::diag((fm::Vector3D){Ixx, Ixx, Izz});
}

Rocket::Rocket(const std::string& filename)
{
	std::fstream fp(filename, std::fstream::in);
	std::string aux_s;
	std::vector<double> aux_v;
	std::vector<std::string> aux_vs;
	uint n_c; // Numero de compomentes
	double Aux;

	// get csv file
	do{
		std::getline(fp, aux_s);
		// while is a comment
	}while(aux_s.find("#") != std::string::npos);
	timeMassThrustDragInput(aux_s);
	
	// get number of components
	do{
		std::getline(fp, aux_s);
		// check if NOT comment
		if(aux_s.find("#") == std::string::npos){
			n_c = std::stoi(aux_s);
			break;
		}
	}while(true);

	Xrb = 0;
	Arb = 0;
	for(uint i = 0; i < n_c; ++i){
		std::getline(fp, aux_s);
		// check if comment
		if(aux_s.find("#") != std::string::npos){
			--i;
			continue;
		}

		// Nose cone moment of inertia
		if(i == 0){
			comp.emplace_back(Component(aux_s, Xrb, cone));
			iTensorBody = comp.back().momInertia();
		} // Initial motor moment of inertia
		else if(i == n_c-1){
			comp.emplace_back(Component(aux_s, Xrb, motor));
			iTensorMotor = comp.back().momInertia();
		} // Other components moment of inertia, cylinder
		else{
			comp.emplace_back(Component(aux_s, Xrb, body));
			iTensorBody += comp.back().momInertia();
		}

		// Rocket dimensions
		Xrb += comp.back().L;
		Aux = M_PI*(fm::sq(comp.back().R));
		Arb = (Aux > Arb) ? Aux : Arb;
	}

	std::cout << "Body Inertia:\n" << iTensorBody << std::endl;
	std::cout << "Motor Initial Inertia:\n" << iTensorMotor << std::endl;
	std::cout << "Total Inertia:\n" << iTensorMotor + iTensorBody << std::endl;

	Xcm = 0;
	Xcm_body = 0;
	Xcm_motor = 0;
	// calculate Xcm of the two groups

	// get fin data
	do{
		std::getline(fp, aux_s);
		// check if NOT comment
		if(aux_s.find("#") == std::string::npos){
			aux_vs = spltString(aux_s);
			nfin = std::stoi(aux_vs[0]);
			ncant = std::stod(aux_vs[1]);
			rf = std::stod(aux_vs[2]);
			Xf = std::stod(aux_vs[3]);
			break;
		}
	}while(true);

	std::cout << nfin << " " << ncant << " " << rf << " " << Xf << "\n";

}

void Rocket::timeMassThrustDragInput(const std::string& filename)
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

fm::Matrix33 Rocket::noseConeInercia(const double& noseConeLength,const double& noseConeRadius, const double& noseConeMass)
{
    double Ixx,Izz;
	
    Ixx = (1.0/10.0) * noseConeMass * noseConeLength * noseConeLength + (3.0/20.0) * noseConeMass * noseConeRadius * noseConeRadius;
    Izz = (3.0/10.0) * noseConeMass * noseConeRadius * noseConeRadius;

	fm::Vector3D inerciaVector = {Ixx,Ixx,Izz};
    return fm::diag(inerciaVector);

}

fm::Matrix33 Rocket::componentInercia(const double& componentLength,const double& componentRadius, const double& componentMass)
{
	double Ixx,Izz;
	
	Ixx = (1.0/12.0) * componentMass * (3 * componentRadius * componentRadius + componentLength * componentLength);
	Izz = (1.0/2.0) * componentMass * (componentRadius * componentRadius);

	fm::Vector3D inerciaVector = {Ixx,Ixx,Izz};
	return fm::diag(inerciaVector);

}











std::vector<std::string> spltString(const std::string& s, const char& cc)
{
    std::vector<std::string> res;
    std::string aux = "";

    for(char c : s){
        if(c == cc){
            if(aux == "") continue;
            res.push_back(aux);
            aux = "";
        }
        else{
            aux += c;
        }
    }
    if(aux != "") res.push_back(aux);

    return res;
}