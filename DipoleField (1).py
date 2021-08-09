import numpy as np
from manim import *
from manim_physics import Charge
import Dipole_Sim


class Field(ThreeDScene):
    def __init__(self, **kwargs):
        self.x_min = -5
        self.x_max = 5
        self.delta_x = 1.5
        self.y_min = -5
        self.y_max = 5
        self.delta_y = 1.5
        self.z_min = -1
        self.z_max = 1
        self.delta_z = 0.7
        self.point_to_field = {}
        self.field = [Dipole_Sim.E, Dipole_Sim.B] #lambda x,t: Dipole_Sim.E(np.array([x[2],1.3,1.7]),0)
        self.field_values = [np.empty((int((self.x_max - self.x_min) // self.delta_x) + 1,
                                     int((self.y_max - self.y_min) // self.delta_y) + 1,
                                     int((self.z_max - self.z_min) // self.delta_z) + 1,
                                     3
                                      )
                                     ) for i in range(2)]
        self.E_COLORS = [TEAL_B, GREEN_B, YELLOW_B, RED_B]
        self.B_COLORS = [BLUE_E, GREEN, YELLOW, PURPLE_D]
        self.correction_factor = 1.75
        self.correction_min_norm = 0.05
        super().__init__(**kwargs)

    def construct(self):
        t = ValueTracker(0)
        # q1 = Charge(Dipole_Sim.q_0, 0.5 * OUT).add_updater(lambda x: x.become(
        #     Charge(
        #         Dipole_Sim.q_0 * np.cos(Dipole_Sim.w * t.get_value()), 0.5 * OUT)))
        # q2 = Charge(-Dipole_Sim.q_0, 0.5 * IN).add_updater(lambda x: x.become(
        #     Charge(
        #         -Dipole_Sim.q_0 * np.cos(Dipole_Sim.w * t.get_value()), 0.5 * OUT)))
        #t_tex = MathTex(f"t={t.get_value()}").to_corner(RIGHT+UP).add_updater(lambda x:x.become(MathTex(f"t={t.get_value()}")))

        self.set_field(self.field[0], t.get_value())
        vg_E = self.get_vg(self.field[0], t.get_value())
        vg_E.add_updater(lambda x: self.update_field(self.field[0], x, t.get_value()))
        self.add(vg_E)
        self.set_camera_orientation(45*DEGREES,45*DEGREES)
        self.set_field(self.field[1], t.get_value())
        vg_B = self.get_vg(self.field[1], t.get_value())
        vg_B.add_updater(lambda x: self.update_field(self.field[1], x, t.get_value()))
        self.add(vg_B)
        #self.wait()

        # self.begin_ambient_camera_rotation()
        self.play(
            t.animate.set_value(0.5),
            run_time=4,
            lag_ratio=0
        )

        #VEC TEST
        # vec = Vector(self.get_field_angle(self.field,np.ones(3),0)[0])
        # self.add(vec,Axes())
        #self.move_camera(90*DEGREES,45*DEGREES, run_time=3)
        # self.play(vec.animate.shift(2*OUT),
        #           vec.animate.rotate(self.get_field_angle(self.field, vec.get_start(), 0)[1],
        #                              get_unit_normal(self.get_field(self.field, vec.get_start(), 0), -OUT),
        #                              vec.get_start()
        #                              )
        #                   )

        # v = VGroup(
        #     *[Vector(Dipole_Sim.B(x,0),).move_to(x) for x in
        #      np.transpose(np.array(np.meshgrid(np.arange(-1,1,0.5),np.arange(-1,1,0.5),np.arange(-1,1,0.5))),
        #                   [1,2,3,0]).reshape((-1,3))
        #     ]
        #     )

    #         self.add(v)
    def get_vg(self, F,t):
        p = 0 if F.__name__ == "E" else 1
        vg = VGroup()
        for z in np.arange(self.z_min, self.z_max, self.delta_z):
            colors = self.E_COLORS if p==0 else self.B_COLORS
            v = ArrowVectorField(lambda x: self.get_field_angle(self.field[p], x, t)[0],
                                 delta_x=self.delta_x,
                                 delta_y=self.delta_y,
                                 x_min=self.x_min,
                                 x_max=self.x_max,
                                 y_min=self.y_min,
                                 y_max=self.y_max,
                                 colors=colors
                                 )
            for vec in v:
                vec.shift(z*OUT)
                vec.rotate(self.get_field_angle(self.field[p], vec.get_start(), t)[1],
                           get_unit_normal(self.get_field(self.field[p], vec.get_start(), 0), OUT),
                           vec.get_start()
                           )
            vg.add(v)
        return vg

    def update_field(self, F, x, t):
        p = 0 if F.__name__ == "E" else 1
        self.set_field(self.field[p], t)
        vg = VGroup()
        for i,z in enumerate(np.arange(self.z_min, self.z_max, self.delta_z)):
            colors = self.E_COLORS if p==0 else self.B_COLORS
            v = ArrowVectorField(lambda x: self.get_field_angle(self.field[p], x, t)[0],
                                 delta_x=self.delta_x,
                                 delta_y=self.delta_y,
                                 x_min=self.x_min,
                                 x_max=self.x_max,
                                 y_min=self.y_min,
                                 y_max=self.y_max,
                                 colors=colors
                                 )
            for j,vec in enumerate(v):
                vec.shift(z*OUT)
                vec.rotate(self.get_field_angle(self.field[p], vec.get_start(), t)[1],
                           get_unit_normal(self.get_field(self.field[p], vec.get_start(), 0), -OUT),
                           vec.get_start()
                           )
                if np.linalg.norm(vec.get_vector())>self.correction_min_norm and np.linalg.norm(vec.get_vector()-x[i][j].get_vector())/np.linalg.norm(vec.get_vector())>self.correction_factor:
                    vec.become(x[i][j])
            vg.add(v)
        x.become(vg)

    def _get_field(self, F, r, t):
        if F.__name__ in self.point_to_field:
            if str((r, t)) in self.point_to_field[F.__name__]:
                return self.point_to_field[F.__name__][str((r, t))]
            else:
                self.point_to_field[F.__name__][str((r, t))] = F(r, t)
                return self.point_to_field[F.__name__][str((r, t))]
        else:
            self.point_to_field[F.__name__] = {}
            self.point_to_field[F.__name__][str((r, t))] = F(r, t)
            return self.point_to_field[F.__name__][str((r, t))]

    def set_field(self, F, t):
        p = 0 if F.__name__=="E" else 1
        for i in np.arange(0, 1 + (self.x_max - self.x_min) // self.delta_x):
            for j in np.arange(0, 1 + (self.y_max - self.y_min) // self.delta_y):
                for k in np.arange(0, 1 + (self.z_max - self.z_min) // self.delta_z):
                    self.field_values[p][int(i)][int(j)][int(k)] = F(
                        np.array([
                            self.x_min + i * self.delta_x,
                            self.y_min + j * self.delta_y,
                            self.z_min + k * self.delta_z
                        ], float),
                        t
                    )

    def get_field(self, F, r, t): #numerical values
        p = 0 if F.__name__ == "E" else 1
        i = (r[0] - self.x_min) // self.delta_x
        j = (r[1] - self.x_min) // self.delta_y
        k = (r[2] - self.z_min) // self.delta_z
        return self.field_values[p][int(i)][int(j)][int(k)]

    def get_field_angle(self, F, r, t):
        f = self.get_field(F, r, t)
        if np.linalg.norm(f) == 0:
            return ORIGIN, 0
        theta = np.pi / 2 - np.arccos(np.dot(f, OUT) / np.linalg.norm(f))
        if np.linalg.norm(np.array([f[0], f[1], 0], float)) == 0:
            return ORIGIN, 0
        vec = (x:=np.array([f[0], f[1], 0], float))/ np.linalg.norm(x) * np.linalg.norm(f)
        # self.point_to_field[r] = theta
        return vec, theta
