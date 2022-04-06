import matplotlib.pyplot as plt
import pandas as pd

dumpfile_cols = ["Particle ID", "Type", "Smoothing length", "Density", "Pressure", "Acceleration", "Velocity", "Position", "Thermal energy"]

def pd_read_dump(filepath):
    df = pd.read_csv(
        filepath, 
        sep = "    ", 
        comment = "#",
        names = dumpfile_cols,
        engine = "python"
    )

    return df

# Density and thermal energy position plots at t = 1
fig, ax = plt.subplots(1, 4, figsize=(20, 5))
df_0 = pd_read_dump("./dumps/0.txt")
df_1 = pd_read_dump("./dumps/200.txt")
print("Read ./dumps/0.txt.")
print("Read ./dumps/200.txt.")

df_alive = df_0[df_0["Type"] == "Alive"]
df_ghost = df_0[df_0["Type"] == "Ghost"]
ax[0].scatter(df_alive["Position"], df_alive["Density"], label="Alive particles", s=60)
ax[0].scatter(df_ghost["Position"], df_ghost["Density"], label="Ghost particles", s=60)
ax[0].set_xlabel("Position")
ax[0].set_ylabel("Density")
ax[0].set_title("Illustration of initial conditions ($T=0$)")

def plot_dump(ax, df, x, y):
    ax.scatter(df[x], df[y], zorder=100)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_xlim(-2, 2)
    ax.set_title(f"{x} vs. {y}")

plot_dump(ax[1], df_1, "Position", "Density")
plot_dump(ax[2], df_1, "Position", "Velocity")
ax[2].set_ylim(-1.5, 1.5)
plot_dump(ax[3], df_1, "Position", "Thermal energy")

fig.legend(loc='upper left')
plt.show()