from __future__ import absolute_import, division, print_function


# helper functions
def get_compartment(data):
    # get data on a single compartment
    skip_keys = [
        'type',
        'simulation_id',
        'experiment_id',
        'time',
        '_id']

    # organize data into a dict
    data_dict = {}
    time_vec = []
    for row in data:
        time = row.get('time')
        time_vec.append(time)
        data_keys = [key for key in row.keys() if key not in skip_keys]
        for key in data_keys:
            value = row[key]
            if key not in data_dict.keys():
                data_dict[key] = {}
            for mol_id, v in value.items():
                if mol_id in data_dict[key].keys():
                    data_dict[key][mol_id].append(v)
                else:
                    data_dict[key][mol_id] = [v]
    data_dict['time'] = time_vec
    return data_dict


def get_lattice(data):
    # get data on the lattice state
    skip_keys = [
        'type',
        'agent_id',
        'simulation_id',
        'experiment_id',
        'time',
        '_id']

    # organize data into a dict
    data_dict = {}
    time_vec = []
    for row in data:
        # get time
        time = row.get('time')
        time_vec.append(time)

        # get agent
        agent_id = row['agent_id']
        if agent_id not in data_dict.keys():
            data_dict[agent_id] = {}

        # get data
        data_keys = [key for key in row.keys() if key not in skip_keys]
        for key in data_keys:
            value = row[key]
            if key in data_dict[agent_id].keys():
                data_dict[agent_id][key].append(value)
            else:
                data_dict[agent_id][key] = [value]

    data_dict['time'] = time_vec
    return data_dict


class Analysis(object):
    def __init__(self, analysis_type='compartment'):
        self.analysis_type = analysis_type

    def requirements(self):
        return []

    def get_data(self, client, query, options={}):
        data = client.find(query)
        return data

    def analyze(self, experiment_config, data):
        pass
