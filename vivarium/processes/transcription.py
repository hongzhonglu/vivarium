from vivarium.actor.process import Process

class Transcription(Process):
    def __init__(self, initial_parameters={}):
        self.promoter_affinities = initial_parameters.get('promoter_affinities', {})

    def next_update(self, timestep, states):
        pass

        
