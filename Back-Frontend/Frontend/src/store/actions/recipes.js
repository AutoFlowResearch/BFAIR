import {
  REQUEST_RECIPES,
  RECEIVE_RECIPES,
  FAILURE_RECIPES,
  ADD_RECIPE_SUCCESS,
  ADD_RECIPE_FAILURE,
  DELETE_RECIPE_SUCCESS,
  UPDATE_RECIPE_SUCCESS
} from "../constants/recipe";
import {
  recipeList,
  createRecipe,
  deleteRecipe as dr,
  updateRecipe as ur
} from "../../services/recipes";

const requestRecipes = () => ({
  type: REQUEST_RECIPES
});

const receiveRecipes = data => ({
  type: RECEIVE_RECIPES,
  payload: data
});

const failureRecipes = error => ({
  type: FAILURE_RECIPES,
  payload: error
});

const addRecipeSuccess = data => ({
  type: ADD_RECIPE_SUCCESS,
  payload: data
});

const addRecipeFailure = error => ({
  type: ADD_RECIPE_FAILURE,
  payload: error
});

export const getAllRecipes = () => async dispatch => {
  dispatch(requestRecipes());
  try {
    const data = await recipeList();
    dispatch(receiveRecipes(data));
  } catch (error) {
    dispatch(failureRecipes(error));
    // throw new Error(error);
  }
};

export const addRecipe = newrecipe => async dispatch => {
  try {
    const data = await createRecipe(newrecipe);
    dispatch(addRecipeSuccess(data));
  } catch (error) {
    console.log("action catch");
    console.log(error.response);
    throw new Error(error);
    // dispatch(addRecipeFailure(error));
  }
};

export const deleteRecipe = id => async dispatch => {
  try {
    const data = await dr(id);
    dispatch({ type: DELETE_RECIPE_SUCCESS, payload: id });
    console.log(data);
  } catch (error) {
    throw new Error(error);
  }
};

export const updateRecipe = (updatedData, id) => async dispatch => {
  try {
    const data = await ur(updatedData, id);
    dispatch({ type: UPDATE_RECIPE_SUCCESS, payload: data });
  } catch (error) {
    throw new Error(error);
  }
};
