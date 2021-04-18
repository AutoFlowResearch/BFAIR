import { combineReducers } from "redux";
import { reducer as formReducer } from "redux-form";
import { auth } from "./auth";
import { usermanagement } from "./users";
import { postManagement } from "./posts";
import { recipeManagement } from "./recipe";

export default combineReducers({
  auth,
  usermanagement,
  postManagement,
  recipeManagement,
  form: formReducer
});
