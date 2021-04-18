import React, { Component } from "react";
import RecipeList from "./RecipeList";
import { Switch } from "react-router-dom";
import { ProtectedRoute } from "../../../common/protectedRoutes";
import NewRecipe from "./NewRecipe";

class Recipe extends Component {
  render() {
    return (
      <div>
        <h1>recipe</h1>
        <RecipeList />
      </div>
    );
  }
}

export default Recipe;
