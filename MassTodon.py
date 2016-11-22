


class MassTodon():
    def __init__( self,
        fastas,
        charges,
        mass_min = 20.,
        mass_max = float('inf')
    ):
        self.fastas     = fastas
        self.charges    = charges
        self.mass_min   = mass_min,
        self.mass_max   = mass_max
