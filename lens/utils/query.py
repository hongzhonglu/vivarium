from __future__ import absolute_import, division, print_function

def get_sims_from_exp(client, experiment_id):
    # given a database client and experiment id, return a list of simulation ids
    simulation_ids = set()
    query_dict = {'experiment_id': experiment_id, 'time': 1.0}
    data = client.find(query_dict)
    for row in data:
        simulation_id = row['simulation_id']
        simulation_ids.add(simulation_id)
    return list(simulation_ids)

def db_to_dict(data):
    # organize data into a dict
    data_dict = {}
    time_vec = []
    for row in data:
        time_vec.append(row['time'])
        data_keys = [key for key in row.keys() if key not in ['simulation_id', 'experiment_id', 'time', '_id']]
        for key in data_keys:
            if key not in data_dict.keys():
                data_dict[key] = {}

            # get data for this key
            key_row = row[key]
            for mol_id, value in key_row.iteritems():


                if mol_id in data_dict[key].keys():
                    data_dict[key][mol_id].append(value)
                else:
                    data_dict[key][mol_id] = [value]
    data_dict['time'] = time_vec

    return data_dict