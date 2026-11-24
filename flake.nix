{
  description = "Protein function prediction: human Swiss-Prot (Track A) + CAFA6 (Track B)";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-26.05";
    flake-parts.url = "github:hercules-ci/flake-parts";
  };

  outputs =
    inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];

      perSystem =
        { pkgs, lib, ... }:
        {
          devShells.default = pkgs.mkShell {
            packages = [
              pkgs.python312
              pkgs.uv
              pkgs.mmseqs2
            ];

            env = lib.optionalAttrs pkgs.stdenv.isLinux {
              # Python wheels load native shared objects with dlopen(3); make the
              # loader aware of the manylinux libraries without RPATH lookups.
              LD_LIBRARY_PATH = lib.makeLibraryPath pkgs.pythonManylinuxPackages.manylinux1;
              # Use the Nix-provided interpreter instead of downloading one.
              UV_PYTHON_DOWNLOADS = "never";
            };

            shellHook = ''
              unset PYTHONPATH
            '';
          };
        };
    };
}
